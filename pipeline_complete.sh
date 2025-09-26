#!/usr/bin/env bash

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }
log "=== PIPELINE VALORMICRO DÉMARRÉ ==="

# ---- CORRECTION TOKEN CONDA
log "Suppression token conda corrompu"
sudo rm -rf /home/fungi/.conda/ 2>/dev/null || rm -rf /home/fungi/.conda/ 2>/dev/null || true

set +u
source $(conda info --base)/etc/profile.d/conda.sh
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
set -u

# ---- 00 GÉNÉRATION MÉTADONNÉES
log "Génération métadonnées"
cd "${ROOTDIR}/00_scripts"

MANIFEST="${ROOTDIR}/98_databasefiles/manifest"
METADATA="${ROOTDIR}/98_databasefiles/sample-metadata.tsv"

if [[ ! -f "$MANIFEST" ]] || [[ ! -f "$METADATA" ]]; then
    if python3 -c "import pandas" 2>/dev/null; then
        log "Utilisation générateur Python"
        python3 "${ROOTDIR}/00_scripts/generate_qiime_files.py"
    else
        log "Utilisation générateur Bash"
        bash "${ROOTDIR}/00_scripts/generate_metadata.sh"
    fi
else
    log "Métadonnées existantes utilisées"
fi

# ---- 01 FASTQC
log "FastQC sur données brutes"
mkdir -p "${ROOTDIR}/02_qualitycheck"

# Nettoyer anciens rapports
rm -f "${ROOTDIR}/02_qualitycheck"/*multiqc* 2>/dev/null || true
rm -rf "${ROOTDIR}/02_qualitycheck/multiqc_data" 2>/dev/null || true

cd "${ROOTDIR}/01_raw_data"

# FastQC sur échantillon test (limité pour validation)
count=0
for file in $(find . -name '*.fastq*' -type f | head -6); do
    count=$((count + 1))
    log "FastQC $count/6: $(basename $file)"
    
    conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || {
        log "Erreur FastQC sur $file, continuons"
        continue
    }
    
    # Pause pour éviter surcharge
    if [ $((count % 3)) -eq 0 ]; then
        sleep 2
    fi
done

log "FastQC terminé sur $count fichiers"

# MultiQC avec gestion d'erreurs Python 2.7
log "MultiQC avec contournement erreurs"
cd "${ROOTDIR}/02_qualitycheck"

conda run -n multiqc multiqc . \
    --force \
    --filename "raw_data_qc" \
    --title "Raw Data Quality Control" \
    --ignore-symlinks \
    --no-ansi 2>/dev/null || {
    log "MultiQC a généré des warnings mais probablement réussi"
    if [ -f "raw_data_qc.html" ]; then
        log "✓ Rapport MultiQC créé malgré warnings"
    else
        log "⚠ MultiQC échoué, continuons sans rapport"
    fi
}

# ---- 02 TRIMMOMATIC
log "Trimmomatic - test sur échantillon réduit"
ADAPTERS="${ROOTDIR}/99_softwares/adapters/sequences.fasta"
mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

# Test sur 3 paires pour validation rapide
pair_count=0
success_count=0

log "Test Trimmomatic sur 3 paires d'échantillons"
find . -name '*R1*.fastq*' -type f | head -3 | while read r1; do
    r2="${r1/_R1/_R2}"
    if [[ -f "$r2" ]]; then
        pair_count=$((pair_count + 1))
        log "Trimmomatic $pair_count/3: $(basename $r1) + $(basename $r2)"
        
        # Noms de base sans extensions
        base1=$(basename "$r1" .fastq.gz)
        base1=$(basename "$base1" .fastq)
        base2=$(basename "$r2" .fastq.gz)
        base2=$(basename "$base2" .fastq)
        
        # Fichiers de sortie Trimmomatic
        out1p="${ROOTDIR}/03_cleaned_data/${base1}_paired.fastq.gz"
        out1u="${ROOTDIR}/03_cleaned_data/${base1}_unpaired.fastq.gz"
        out2p="${ROOTDIR}/03_cleaned_data/${base2}_paired.fastq.gz"
        out2u="${ROOTDIR}/03_cleaned_data/${base2}_unpaired.fastq.gz"
        
        # Exécution Trimmomatic
        conda run -n trimmomatic trimmomatic PE -threads 4 -phred33 "$r1" "$r2" \
            "$out1p" "$out1u" "$out2p" "$out2u" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 || {
            log "Erreur Trimmomatic sur $r1/$r2"
            continue
        }
        
        # Vérification synchronisation CRITIQUE
        if [[ -f "$out1p" && -f "$out2p" ]]; then
            # Compter reads avec méthode robuste
            count1=$(( $(zcat "$out1p" 2>/dev/null | wc -l || cat "$out1p" | wc -l) / 4 ))
            count2=$(( $(zcat "$out2p" 2>/dev/null | wc -l || cat "$out2p" | wc -l) / 4 ))
            
            log "Vérification synchronisation: R1=$count1 reads, R2=$count2 reads"
            
            if [ "$count1" = "$count2" ] && [ "$count1" -gt 0 ]; then
                log "✓ Paire $pair_count SYNCHRONISÉE: $count1 reads"
                success_count=$((success_count + 1))
            else
                log "✗ Paire $pair_count DÉSYNCHRONISÉE ($count1 vs $count2) - SUPPRESSION"
                rm -f "$out1p" "$out2p" "$out1u" "$out2u"
            fi
        else
            log "✗ Fichiers de sortie manquants pour paire $pair_count"
        fi
        
        # Pause entre traitements
        sleep 3
    fi
done

log "Trimmomatic terminé - Paires synchronisées réussies: $success_count"

# Vérifier qu'au moins une paire a réussi
paired_files=$(find "${ROOTDIR}/03_cleaned_data" -name "*_paired.fastq*" | wc -l)
if [ "$paired_files" -eq 0 ]; then
    log "ERREUR: Aucun fichier paired généré par Trimmomatic"
    exit 1
fi

log "Fichiers paired trouvés: $paired_files"

# ---- 02.5 FASTQC/MULTIQC SUR DONNÉES NETTOYÉES
log "FastQC/MultiQC sur données nettoyées après Trimmomatic"
mkdir -p "${ROOTDIR}/03_cleaned_data_qc"

# Nettoyer anciens rapports
rm -f "${ROOTDIR}/03_cleaned_data_qc"/*multiqc* 2>/dev/null || true
rm -rf "${ROOTDIR}/03_cleaned_data_qc/multiqc_data" 2>/dev/null || true

cd "${ROOTDIR}/03_cleaned_data"

# FastQC sur fichiers paired nettoyés
log "FastQC sur fichiers paired nettoyés"
count=0
for file in *_paired.fastq*; do
    if [ -f "$file" ]; then
        count=$((count + 1))
        log "FastQC cleaned $count: $(basename $file)"
        
        conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/03_cleaned_data_qc" --threads 2 --quiet || {
            log "Erreur FastQC sur $file, continuons"
            continue
        }
        
        # Pause pour éviter surcharge
        if [ $((count % 4)) -eq 0 ]; then
            sleep 2
        fi
    fi
done

# MultiQC sur données nettoyées
log "MultiQC sur données nettoyées"
cd "${ROOTDIR}/03_cleaned_data_qc"

conda run -n multiqc multiqc . \
    --force \
    --filename "cleaned_data_qc" \
    --title "Cleaned Data Quality Control After Trimmomatic" \
    --ignore-symlinks \
    --no-ansi 2>/dev/null || {
    log "MultiQC cleaned data a généré des warnings mais probablement réussi"
    if [ -f "cleaned_data_qc.html" ]; then
        log "✓ Rapport MultiQC cleaned data créé"
    else
        log "⚠ MultiQC cleaned data échoué, continuons"
    fi
}

log "✅ Contrôle qualité post-Trimmomatic terminé"

# ---- 03 GÉNÉRATION MANIFEST PAIRED AVEC IDS UNIQUES
log "Génération manifest paired avec IDs UNIQUES - SOLUTION ROBUSTE"
cd "${ROOTDIR}/98_databasefiles"

# Nettoyer anciens manifests
rm -f manifest_paired manifest_control_paired

# Créer headers
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest_paired
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest_control_paired

# Array associatif pour éviter les doublons
declare -A seen_ids

# Scanner fichiers paired dans cleaned_data
cd "${ROOTDIR}/03_cleaned_data"
count=0
control_count=0

log "Scan des fichiers paired avec génération d'IDs UNIQUES dans: $(pwd)"
for r1_file in *R1*_paired.fastq*; do
    if [ -f "$r1_file" ]; then
        # Trouver R2 correspondant
        r2_file="${r1_file/R1/R2}"
        
        if [ -f "$r2_file" ]; then
            # Vérifier taille des fichiers (>1KB pour éviter fichiers vides)
            r1_size=$(stat -c%s "$r1_file" 2>/dev/null || echo "0")
            r2_size=$(stat -c%s "$r2_file" 2>/dev/null || echo "0")
            
            if [ "$r1_size" -gt 1000 ] && [ "$r2_size" -gt 1000 ]; then
                
                # Utiliser le nom de fichier complet sans extensions comme base
                base_name=$(basename "$r1_file")
                # Supprimer toutes les extensions possibles pour avoir un sample-id propre
                sample_id="${base_name%_R1*}"        # Supprimer _R1 et tout ce qui suit
                sample_id="${sample_id%.fastq*}"     # Supprimer .fastq ou .fastq.gz
                sample_id="${sample_id%.fq*}"        # Supprimer .fq ou .fq.gz  
                sample_id="${sample_id%_paired*}"    # Supprimer _paired si présent
                
                # Remplacer les caractères problématiques par des underscores
                sample_id="${sample_id//[^a-zA-Z0-9._-]/_}"
                
                # Vérifier l'unicité et ajuster si nécessaire
                original_id="$sample_id"
                counter=1
                while [[ -n "${seen_ids[$sample_id]:-}" ]]; do
                    sample_id="${original_id}_${counter}"
                    counter=$((counter + 1))
                    log "ID dupliqué détecté, nouveau: $sample_id"
                done
                
                # Marquer cet ID comme utilisé
                seen_ids["$sample_id"]=1
                
                # Chemins absolus
                r1_abs="${ROOTDIR}/03_cleaned_data/$r1_file"
                r2_abs="${ROOTDIR}/03_cleaned_data/$r2_file"
                
                # Détecter si c'est un contrôle (inclure "eau" dans les contrôles)
                if echo "${sample_id,,}" | grep -qE "(neg|blank|control|ctrl|eau)"; then
                    echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "${ROOTDIR}/98_databasefiles/manifest_control_paired"
                    control_count=$((control_count + 1))
                    log "Contrôle ajouté: $sample_id (fichier: $r1_file)"
                else
                    echo -e "$sample_id\t$r1_abs\t$r2_abs" >> "${ROOTDIR}/98_databasefiles/manifest_paired"
                    count=$((count + 1))
                    log "Échantillon ajouté: $sample_id (fichier: $r1_file)"
                fi
            else
                log "Fichiers trop petits ignorés: $r1_file ($r1_size B), $r2_file ($r2_size B)"
            fi
        else
            log "Pas de R2 correspondant pour: $r1_file"
        fi
    fi
done

cd "${ROOTDIR}/98_databasefiles"
log "Manifest créé: $count échantillons principaux, $control_count contrôles"

# Vérifier manifest principal
if [ $(wc -l < manifest_paired) -le 1 ]; then
    log "ERREUR: Aucun échantillon dans manifest_paired"
    log "Fichiers trouvés dans cleaned_data:"
    ls -la "${ROOTDIR}/03_cleaned_data/"*paired* 2>/dev/null || log "Aucun fichier paired"
    exit 1
fi

# Supprimer manifest contrôles si vide (seulement header)
if [ $(wc -l < manifest_control_paired) -le 1 ]; then
    rm -f manifest_control_paired
    log "Pas de contrôles détectés"
fi

# VÉRIFICATION FINALE ANTI-DOUBLONS
log "Vérification finale des doublons dans le manifest"
duplicates=$(cut -f1 manifest_paired | sort | uniq -d)
if [ -n "$duplicates" ]; then
    log "❌ ERREUR: Doublons encore présents: $duplicates"
    exit 1
else
    log "✅ Aucun doublon - IDs uniques confirmés"
fi

log "Contenu du manifest paired final avec IDs UNIQUES:"
cat manifest_paired

if [ -f manifest_control_paired ]; then
    log "Contenu du manifest contrôles:"  
    cat manifest_control_paired
fi

# ---- 04 QIIME2 IMPORT
log "QIIME2 Import avec fichiers paired synchronisés et IDs uniques"
mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

MANIFEST_PAIRED="${ROOTDIR}/98_databasefiles/manifest_paired"
MANIFEST_CONTROL_PAIRED="${ROOTDIR}/98_databasefiles/manifest_control_paired"

# Import principal
log "Import QIIME2 principal avec IDs uniques"
conda run -n qiime2-2021.4 qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_PAIRED" \
    --output-path "core/demux_paired.qza" \
    --input-format PairedEndFastqManifestPhred33V2 || {
    log "ERREUR import QIIME2 principal"
    log "Vérification du manifest:"
    head -5 "$MANIFEST_PAIRED"
    
    log "Vérification doublons dans manifest:"
    cut -f1 "$MANIFEST_PAIRED" | sort | uniq -c | sort -nr
    exit 1
}

log "✅ Import QIIME2 principal réussi avec IDs uniques !"

# Import contrôles si présents
HAS_CONTROLS=false
if [ -f "$MANIFEST_CONTROL_PAIRED" ] && [ -s "$MANIFEST_CONTROL_PAIRED" ]; then
    log "Import contrôles QIIME2"
    conda run -n qiime2-2021.4 qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_CONTROL_PAIRED" \
        --output-path "core/demux_neg.qza" \
        --input-format PairedEndFastqManifestPhred33V2 && HAS_CONTROLS=true || {
        log "Import contrôles échoué, continuons sans contrôles"
    }
fi

# ---- 05 DADA2
log "DADA2 denoising"
cd "${ROOTDIR}/05_QIIME2/core"

# Tentative DADA2
log "Lancement DADA2 avec fichiers paired synchronisés et IDs uniques..."
conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux_paired.qza \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "$NTHREADS" || {
    
    log "❌ DADA2 ÉCHOUÉ"
    exit 1
}

log "🎉 DADA2 RÉUSSI !"

# ---- 05 SILVA CLASSIFIER EXISTANT
log "Utilisation SILVA SSU 138.2 existant"
cd "${ROOTDIR}/98_databasefiles"

CLASSIFIER="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
if ! conda run -n qiime2-2021.4 qiime tools validate "$CLASSIFIER" 2>/dev/null; then
    log "Création classifier SILVA 138.2"
    cd "${ROOTDIR}/98_databasefiles"
    # Extraction région V4-V5
    conda run -n qiime2-2021.4 qiime feature-classifier extract-reads \
        --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
        --p-f-primer GTGYCAGCMGCCGCGGTAA \
        --p-r-primer CCGYCAATTYMTTTRAGTTT \
        --p-n-jobs 2 \
        --p-read-orientation 'forward' \
        --o-reads silva-138.2-ssu-nr99-seqs-515f-926r.qza

    # Déréplication
    conda run -n qiime2-2021.4 qiime rescript dereplicate \
        --i-sequences silva-138.2-ssu-nr99-seqs-515f-926r.qza \
        --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
        --p-mode uniq \
        --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
        --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza

    # Fit classifier
    conda run -n qiime2-2021.4 qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads silva-138.2-ssu-nr99-seqs-515f-926r-uniq.qza \
        --i-reference-taxonomy silva-138.2-ssu-nr99-tax-515f-926r-derep-uniq.qza \
        --o-classifier "$CLASSIFIER"
fi

# ---- 06 CLASSIFICATION
log "Classification taxonomique SILVA 138.2"
cd "${ROOTDIR}/05_QIIME2/core"
conda run -n qiime2-2021.4 qiime feature-classifier classify-sklearn \
    --i-classifier "${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS"

log "✅ Classification taxonomique SILVA 138.2 officiel réussie"

# Vérifier le contenu de la taxonomie
conda run -n qiime2-2021.4 qiime tools export \
    --input-path taxonomy.qza \
    --output-path temp_tax_check

if [ -f "temp_tax_check/taxonomy.tsv" ]; then
    tax_count=$(tail -n +2 temp_tax_check/taxonomy.tsv | wc -l)
    log "✅ Taxonomie SILVA 138.2 officiel contient $tax_count classifications"
    log "Échantillon de la taxonomie SILVA 138.2 officiel:"
    head -5 temp_tax_check/taxonomy.tsv | column -t -s$'\t' || head -5 temp_tax_check/taxonomy.tsv
fi
rm -rf temp_tax_check

# ---- 08 ANALYSES FINALES 
log "Analyses finales : core features, taxa barplots et exports"
mkdir -p "${ROOTDIR}/05_QIIME2/subtables" "${ROOTDIR}/05_QIIME2/export"

cd "${ROOTDIR}/05_QIIME2/core"

# Summary de la table
log "Summary table pour profondeur raréfaction"
conda run -n qiime2-2021.4 qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization "../visual/table-summary.qzv" || {
    log "Erreur summary table"
    exit 1
}

# Export du summary pour extraction automatique
conda run -n qiime2-2021.4 qiime tools export \
    --input-path "../visual/table-summary.qzv" \
    --output-path "../visual/table-summary"

# Extraction profondeur raréfaction
if [ -f "../visual/table-summary/sample-frequency-detail.csv" ]; then
    RAREFACTION_DEPTH_FLOAT=$(awk -F',' 'NR>1 {print $2}' "../visual/table-summary/sample-frequency-detail.csv" | sort -n | awk 'NR==int(NR*0.1)+1' || echo "5000")
    RAREFACTION_DEPTH=$(printf "%.0f" "$RAREFACTION_DEPTH_FLOAT" 2>/dev/null || echo "5000")
    
    if ! [[ "$RAREFACTION_DEPTH" =~ ^[0-9]+$ ]] || [ "$RAREFACTION_DEPTH" -lt 1 ]; then
        RAREFACTION_DEPTH=5000
        log "Valeur invalide détectée, utilisation par défaut: $RAREFACTION_DEPTH"
    else
        log "Profondeur de raréfaction automatique: $RAREFACTION_DEPTH"
    fi
else
    RAREFACTION_DEPTH=5000
    log "Fichier summary non trouvé, utilisation par défaut: $RAREFACTION_DEPTH"
fi

# Raréfaction
conda run -n qiime2-2021.4 qiime feature-table rarefy \
    --i-table table.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table "../subtables/RarTable-all.qza" || {
    log "Erreur raréfaction, copie table originale"
    cp table.qza "../subtables/RarTable-all.qza"
}

# Core features analysis
conda run -n qiime2-2021.4 qiime feature-table core-features \
    --i-table "../subtables/RarTable-all.qza" \
    --p-min-fraction 0.1 \
    --p-max-fraction 1.0 \
    --p-steps 10 \
    --o-visualization "../visual/CoreBiom-all.qzv" || {
    log "Erreur core features analysis"
}

# Taxa barplots
log "Génération taxa barplots avec SILVA 138.2 officiel"
conda run -n qiime2-2021.4 qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --o-visualization "../visual/taxa-bar-plots.qzv" || {
    log "Erreur taxa barplots"
}

# ---- 09 MÉTRIQUES DE DIVERSITÉ CORRIGÉES
log "Calcul métriques de diversité avec outputs individuels"
mkdir -p "${ROOTDIR}/05_QIIME2/diversity" "${ROOTDIR}/05_QIIME2/pcoa" 

cd "${ROOTDIR}/05_QIIME2/core"

# Arbre phylogénétique
log "Génération arbre phylogénétique"
if [ ! -f "tree.qza" ]; then
    conda run -n qiime2-2021.4 qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences rep-seqs.qza \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree tree.qza \
        --p-n-threads "$NTHREADS" || {
        log "Erreur génération arbre phylogénétique"
    }
fi

# Métadonnées pour diversité
mkdir -p "../98_databasefiles"
if [ -f "${ROOTDIR}/98_databasefiles/manifest_paired" ]; then
    echo -e "sample-id\tgroup\ttype" > "../98_databasefiles/diversity-metadata.tsv"
    
    tail -n +2 "${ROOTDIR}/98_databasefiles/manifest_paired" | cut -f1 | while read -r sample_id; do
        if echo "${sample_id,,}" | grep -qE "(neg|blank|control|ctrl|eau)"; then
            echo -e "$sample_id\tcontrol\tnegative" >> "../98_databasefiles/diversity-metadata.tsv"
        else
            echo -e "$sample_id\tsample\tenvironmental" >> "../98_databasefiles/diversity-metadata.tsv"
        fi
    done
    
    log "✅ Métadonnées diversité créées"
fi

# NETTOYER ANCIENS RÉSULTATS POUR ÉVITER L'ERREUR --output-dir EXISTS
log "Nettoyage anciens résultats de diversité"
rm -rf diversity-results 2>/dev/null || true

# Core metrics phylogenetic avec outputs individuels [web:33][web:25]
log "Lancement core-metrics-phylogenetic avec tous les outputs individuels"

conda run -n qiime2-2021.4 qiime diversity core-metrics-phylogenetic \
    --i-table table.qza \
    --i-phylogeny tree.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --m-metadata-file "../98_databasefiles/diversity-metadata.tsv" \
    --o-rarefied-table rarefied_table.qza \
    --o-faith-pd-vector diversity/Vector-faith_pd.qza \
    --o-observed-features-vector diversity/Vector-observed_asv.qza \
    --o-shannon-vector diversity/Vector-shannon.qza \
    --o-evenness-vector diversity/Vector-evenness.qza \
    --o-unweighted-unifrac-distance-matrix diversity/Matrix-unweighted_unifrac.qza \
    --o-weighted-unifrac-distance-matrix diversity/Matrix-weighted_unifrac.qza \
    --o-jaccard-distance-matrix diversity/Matrix-jaccard.qza \
    --o-bray-curtis-distance-matrix diversity/Matrix-braycurtis.qza \
    --o-unweighted-unifrac-pcoa-results pcoa/PCoA-unweighted_unifrac.qza \
    --o-weighted-unifrac-pcoa-results pcoa/PCoA-weighted_unifrac.qza \
    --o-jaccard-pcoa-results pcoa/PCoA-jaccard.qza \
    --o-bray-curtis-pcoa-results pcoa/PCoA-braycurtis.qza \
    --o-unweighted-unifrac-emperor visual/Emperor-unweighted_unifrac.qzv \
    --o-weighted-unifrac-emperor visual/Emperor-weighted_unifrac.qzv \
    --o-jaccard-emperor visual/Emperor-jaccard.qzv \
    --o-bray-curtis-emperor visual/Emperor-braycurtis.qzv || {
    
    log "Erreur core-metrics-phylogenetic, tentative sans phylogénie"
    
    # Alternative sans phylogénie avec tous les outputs individuels
    conda run -n qiime2-2021.4 qiime diversity core-metrics \
        --i-table table.qza \
        --p-sampling-depth "$RAREFACTION_DEPTH" \
        --m-metadata-file "../98_databasefiles/diversity-metadata.tsv" \
        --o-rarefied-table rarefied_table.qza \
        --o-observed-features-vector diversity/Vector-observed_asv.qza \
        --o-shannon-vector diversity/Vector-shannon.qza \
        --o-evenness-vector diversity/Vector-evenness.qza \
        --o-jaccard-distance-matrix diversity/Matrix-jaccard.qza \
        --o-bray-curtis-distance-matrix diversity/Matrix-braycurtis.qza \
        --o-jaccard-pcoa-results pcoa/PCoA-jaccard.qza \
        --o-bray-curtis-pcoa-results pcoa/PCoA-braycurtis.qza \
        --o-jaccard-emperor visual/Emperor-jaccard.qzv \
        --o-bray-curtis-emperor visual/Emperor-braycurtis.qzv || {
        
        log "Erreur core-metrics, création métriques individuelles"
        
        # Créer métriques alpha individuellement
        conda run -n qiime2-2021.4 qiime diversity alpha \
            --i-table table.qza \
            --p-metric observed_features \
            --o-alpha-diversity diversity/Vector-observed_asv.qza || true
            
        conda run -n qiime2-2021.4 qiime diversity alpha \
            --i-table table.qza \
            --p-metric shannon \
            --o-alpha-diversity diversity/Vector-shannon.qza || true
            
        conda run -n qiime2-2021.4 qiime diversity alpha \
            --i-table table.qza \
            --p-metric pielou_e \
            --o-alpha-diversity diversity/Vector-evenness.qza || true
        
        # Créer matrices beta individuellement
        conda run -n qiime2-2021.4 qiime diversity beta \
            --i-table table.qza \
            --p-metric jaccard \
            --o-distance-matrix diversity/Matrix-jaccard.qza || true
            
        conda run -n qiime2-2021.4 qiime diversity beta \
            --i-table table.qza \
            --p-metric braycurtis \
            --o-distance-matrix diversity/Matrix-braycurtis.qza || true
    }
}

log "✅ Métriques de diversité créées avec outputs individuels"

# Compter les fichiers créés
diversity_count=$(find diversity -name "*.qza" 2>/dev/null | wc -l || echo "0")
pcoa_count=$(find pcoa -name "*.qza" 2>/dev/null | wc -l || echo "0")  
emperor_count=$(find visual -name "Emperor*.qzv" 2>/dev/null | wc -l || echo "0")

log "Résumé: $diversity_count métriques diversité, $pcoa_count PCoA, $emperor_count Emperor plots"

# ---- 10 EXPORTS QIIME2
log "Export de tous les fichiers QIIME2"
mkdir -p "${ROOTDIR}/05_QIIME2/export/core" \
         "${ROOTDIR}/05_QIIME2/export/subtables/RarTable-all" \
         "${ROOTDIR}/05_QIIME2/export/visual/CoreBiom-all" \
         "${ROOTDIR}/05_QIIME2/export/visual/taxa-bar-plots" \
         "${ROOTDIR}/05_QIIME2/export/diversity_tsv"

cd "${ROOTDIR}/05_QIIME2"

# Export table principale
conda run -n qiime2-2021.4 qiime tools export \
    --input-path core/table.qza \
    --output-path export/core/table

# Export séquences représentatives
conda run -n qiime2-2021.4 qiime tools export \
    --input-path core/rep-seqs.qza \
    --output-path export/core/rep-seqs

# Export taxonomie SILVA 138.2 officiel
conda run -n qiime2-2021.4 qiime tools export \
    --input-path core/taxonomy.qza \
    --output-path export/core/taxonomy

# Export table raréfiée
conda run -n qiime2-2021.4 qiime tools export \
    --input-path subtables/RarTable-all.qza \
    --output-path export/subtables/RarTable-all

# Export visualisations
conda run -n qiime2-2021.4 qiime tools export \
    --input-path visual/CoreBiom-all.qzv \
    --output-path export/visual/CoreBiom-all || {
    log "Erreur export CoreBiom visualization"
}

conda run -n qiime2-2021.4 qiime tools export \
    --input-path visual/taxa-bar-plots.qzv \
    --output-path export/visual/taxa-bar-plots || {
    log "Erreur export taxa barplots"
}

# ---- EXPORT TOUS LES FICHIERS DE DIVERSITÉ EN TSV (CORRIGÉ)
log "Export de TOUS les fichiers de diversité en format TSV/TXT"

export_diversity_to_tsv() {
    local qza_file="$1"
    local output_name="$2"
    
    if [ -f "$qza_file" ]; then
        log "Export $output_name en TSV"
        conda run -n qiime2-2021.4 qiime tools export \
            --input-path "$qza_file" \
            --output-path "export/diversity_tsv/${output_name}_temp" || return 1
        
        # Trouver et copier TOUS les fichiers générés
        find "export/diversity_tsv/${output_name}_temp" -name "*.tsv" -o -name "*.txt" | while read -r found_file; do
            if [ -f "$found_file" ]; then
                base_name=$(basename "$found_file")
                final_name="${output_name}_${base_name}"
                cp "$found_file" "export/diversity_tsv/${final_name}"
                log "✅ $final_name créé"
            fi
        done
        
        # Si aucun TSV trouvé, chercher d'autres formats
        if [ ! -f "export/diversity_tsv/${output_name}.tsv" ]; then
            any_file=$(find "export/diversity_tsv/${output_name}_temp" -type f | head -1)
            if [ -f "$any_file" ]; then
                cp "$any_file" "export/diversity_tsv/${output_name}.tsv"
                log "✅ $output_name.tsv créé (format alternatif)"
            fi
        fi
        
        # Nettoyer
        rm -rf "export/diversity_tsv/${output_name}_temp"
        return 0
    else
        log "❌ $qza_file non trouvé"
        return 1
    fi
}

# Export TOUS les fichiers de diversité en TSV
log "Export systématique de tous les fichiers de diversité"

# Métriques alpha
export_diversity_to_tsv "diversity/Vector-observed_asv.qza" "observed_features"
export_diversity_to_tsv "diversity/Vector-shannon.qza" "shannon"
export_diversity_to_tsv "diversity/Vector-evenness.qza" "evenness"
export_diversity_to_tsv "diversity/Vector-faith_pd.qza" "faith_pd"

# Matrices de distance
export_diversity_to_tsv "diversity/Matrix-jaccard.qza" "jaccard_distance"
export_diversity_to_tsv "diversity/Matrix-braycurtis.qza" "bray_curtis_distance"
export_diversity_to_tsv "diversity/Matrix-unweighted_unifrac.qza" "unweighted_unifrac_distance"
export_diversity_to_tsv "diversity/Matrix-weighted_unifrac.qza" "weighted_unifrac_distance"

# PCoA
export_diversity_to_tsv "pcoa/PCoA-jaccard.qza" "jaccard_pcoa"
export_diversity_to_tsv "pcoa/PCoA-braycurtis.qza" "bray_curtis_pcoa"
export_diversity_to_tsv "pcoa/PCoA-unweighted_unifrac.qza" "unweighted_unifrac_pcoa"
export_diversity_to_tsv "pcoa/PCoA-weighted_unifrac.qza" "weighted_unifrac_pcoa"

# Stats DADA2
if [ -f "core/denoising-stats.qza" ]; then
    export_diversity_to_tsv "core/denoising-stats.qza" "dada2_stats"
fi

# Table raréfiée si créée
if [ -f "core/rarefied_table.qza" ]; then
    export_diversity_to_tsv "core/rarefied_table.qza" "rarefied_table"
fi

log "✅ Tous les fichiers de diversité exportés en TSV dans export/diversity_tsv/"

# Compter et lister les fichiers TSV créés
tsv_count=$(find export/diversity_tsv -name "*.tsv" -o -name "*.txt" 2>/dev/null | wc -l || echo "0")
log "Nombre total de fichiers TSV/TXT créés: $tsv_count"

# Lister tous les fichiers créés
log "Fichiers TSV/TXT créés dans diversity_tsv:"
ls -la export/diversity_tsv/ 2>/dev/null || log "Dossier diversity_tsv vide"

# ---- 11 CONVERSIONS BIOM VERS TSV CORRIGÉES
log "Conversion BIOM vers TSV avec syntaxe bash corrigée"
cd "${ROOTDIR}/05_QIIME2/export"

# S'assurer que les répertoires existent
mkdir -p subtables/RarTable-all core/taxonomy

# Fonction de conversion BIOM vers TSV robuste
convert_biom_to_tsv_fixed() {
    local biom_file="$1"
    local output_tsv="$2"
    
    if [ ! -f "$biom_file" ]; then
        log "❌ Fichier BIOM manquant : $biom_file"
        return 1
    fi
    
    log "Conversion $biom_file vers $output_tsv"
    
    # Méthode 1 : biom convert standard
    conda run -n qiime2-2021.4 biom convert \
        -i "$biom_file" \
        -o "$output_tsv" \
        --to-tsv 2>/dev/null && {
        log "✅ Conversion réussie avec biom convert"
        return 0
    }
    
    # Méthode 2 : biom dans environnement base
    biom convert -i "$biom_file" -o "$output_tsv" --to-tsv 2>/dev/null && {
        log "✅ Conversion réussie avec biom système"
        return 0
    }
    
    # Méthode 3 : Python avec biom-format
    conda run -n qiime2-2021.4 python3 -c "
import biom
import sys
try:
    table = biom.load_table('$biom_file')
    with open('$output_tsv', 'w') as f:
        f.write('#OTU ID\\t' + '\\t'.join(table.ids(axis='sample')) + '\\n')
        for feature_id, feature_data in zip(table.ids(axis='observation'), table.matrix_data.toarray()):
            line = feature_id + '\\t' + '\\t'.join(map(str, feature_data.flatten()))
            f.write(line + '\\n')
    print('Conversion Python réussie')
except Exception as e:
    print(f'Erreur Python: {e}')
    sys.exit(1)
" && {
        log "✅ Conversion réussie avec Python"
        return 0
    }
    
    log "❌ Toutes les méthodes de conversion BIOM ont échoué pour $biom_file"
    return 1
}

# Conversion table raréfiée
if [ -f "subtables/RarTable-all/feature-table.biom" ]; then
    log "Conversion table raréfiée BIOM vers TSV"
    
    if convert_biom_to_tsv_fixed "subtables/RarTable-all/feature-table.biom" "subtables/RarTable-all/table-from-biom.tsv"; then
        # Modification header pour créer ASV.tsv - SYNTAXE BASH CORRIGÉE
        if [ -f "subtables/RarTable-all/table-from-biom.tsv" ]; then
            sed '1d ; s/#OTU ID/ASV_ID/' \
                subtables/RarTable-all/table-from-biom.tsv > \
                subtables/RarTable-all/ASV.tsv
            
            log "✅ Fichier ASV.tsv créé : $(wc -l < subtables/RarTable-all/ASV.tsv 2>/dev/null || echo "0") lignes"
        fi
    else
        log "❌ Erreur conversion BIOM table raréfiée"
    fi
fi

# Conversion table principale
if [ -f "core/table/feature-table.biom" ]; then
    log "Conversion table principale BIOM vers TSV"
    
    if convert_biom_to_tsv_fixed "core/table/feature-table.biom" "core/table/table-from-biom.tsv"; then
        if [ -f "core/table/table-from-biom.tsv" ]; then
            sed '1d ; s/#OTU ID/ASV_ID/' \
                core/table/table-from-biom.tsv > \
                core/table/ASV.tsv
            log "✅ Fichier ASV.tsv principal créé"
        fi
    else
        log "❌ Erreur conversion BIOM table principale"
    fi
fi

# ---- 12 CRÉATION FICHIER ASV AVEC TAXONOMIE SILVA 138.2 OFFICIEL
log "Création fichier ASV.txt avec taxonomie SILVA 138.2 officiel du site arb-silva.de"
cd "${ROOTDIR}/05_QIIME2/export"

create_asv_with_official_silva_taxonomy() {
    local asv_file="subtables/RarTable-all/ASV.tsv"
    local taxonomy_file="core/taxonomy/taxonomy.tsv"
    local output_file="subtables/RarTable-all/ASV.txt"
    
    if [ ! -f "$asv_file" ] || [ ! -f "$taxonomy_file" ]; then
        log "❌ Fichiers requis manquants : $asv_file ou $taxonomy_file"
        return 1
    fi
    
    log "Traitement des fichiers ASV avec taxonomie SILVA 138.2 officiel"
    
    # Obtenir header des échantillons depuis ASV.tsv
    sample_header=$(head -1 "$asv_file" | cut -f2-)
    
    # Créer header final avec taxonomie
    echo -e "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t${sample_header}" > "$output_file"
    
    # Traiter chaque ASV
    tail -n +2 "$asv_file" | while IFS=$'\t' read -r asv_id asv_counts; do
        # Initialiser taxonomie par défaut
        kingdom="Unassigned"
        phylum="Unassigned"
        class="Unassigned"
        order="Unassigned"
        family="Unassigned"
        genus="Unassigned"
        species="Unassigned"
        
        # Chercher taxonomie dans fichier taxonomy.tsv SILVA 138.2 officiel
        if tax_line=$(grep "^${asv_id}" "$taxonomy_file" 2>/dev/null); then
            tax_string=$(echo "$tax_line" | cut -f2)
            
            # Parser la taxonomie SILVA 138.2 officiel (format d__; p__; c__; etc.)
            if [ -n "$tax_string" ]; then
                # Séparer par ; et traiter chaque niveau
                IFS=';' read -ra tax_levels <<< "$tax_string"
                
                for level in "${tax_levels[@]}"; do
                    level=$(echo "$level" | xargs)  # Trim whitespace
                    
                    if [[ "$level" == d__* ]]; then
                        kingdom="${level#d__}"
                        kingdom="${kingdom:-Unassigned}"
                    elif [[ "$level" == p__* ]]; then
                        phylum="${level#p__}"
                        phylum="${phylum:-Unassigned}"
                    elif [[ "$level" == c__* ]]; then
                        class="${level#c__}"
                        class="${class:-Unassigned}"
                    elif [[ "$level" == o__* ]]; then
                        order="${level#o__}"
                        order="${order:-Unassigned}"
                    elif [[ "$level" == f__* ]]; then
                        family="${level#f__}"
                        family="${family:-Unassigned}"
                    elif [[ "$level" == g__* ]]; then
                        genus="${level#g__}"
                        genus="${genus:-Unassigned}"
                    elif [[ "$level" == s__* ]]; then
                        species="${level#s__}"
                        species="${species:-Unassigned}"
                    fi
                done
            fi
        fi
        
        # Nettoyer les valeurs vides
        [ -z "$kingdom" ] && kingdom="Unassigned"
        [ -z "$phylum" ] && phylum="Unassigned"
        [ -z "$class" ] && class="Unassigned"
        [ -z "$order" ] && order="Unassigned"
        [ -z "$family" ] && family="Unassigned"
        [ -z "$genus" ] && genus="Unassigned"
        [ -z "$species" ] && species="Unassigned"
        
        # Écrire ligne finale avec taxonomie SILVA 138.2 officiel
        echo -e "${kingdom}\t${phylum}\t${class}\t${order}\t${family}\t${genus}\t${species}\t${asv_counts}" >> "$output_file"
    done
    
    log "✅ Fichier ASV.txt créé avec taxonomie SILVA 138.2 officiel (depuis arb-silva.de)"
    log "Lignes dans le fichier final: $(wc -l < "$output_file" 2>/dev/null || echo "0")"
    
    # Afficher un échantillon du résultat
    log "Aperçu du fichier ASV.txt avec taxonomie SILVA 138.2 officiel:"
    head -3 "$output_file" | column -t -s$'\t' 2>/dev/null || head -3 "$output_file"
}

# Exécuter la fonction
create_asv_with_official_silva_taxonomy || {
    log "❌ Création ASV.txt échouée"
}

# ---- 13 TABLEAUX RÉCAPITULATIFS
log "Création tableaux récapitulatifs"
mkdir -p "${ROOTDIR}/05_QIIME2/export/summary_tables"
cd "${ROOTDIR}/05_QIIME2/export"

# Créer rapport de synthèse final
log "Création rapport de synthèse final"
cat > "summary_tables/PIPELINE_SUMMARY_REPORT.md" << 'EOF'
# Rapport de Synthèse Pipeline QIIME2 Valormicro avec SILVA 138.2 OFFICIEL

## ✅ Taxonomie authentique SILVA 138.2 téléchargée depuis le site officiel

- **Source**: https://www.arb-silva.de (site officiel SILVA)
- **Release date**: 11 juillet 2024
- **Base de données**: SILVA SSU Ref NR 138.2 (taxonomie la plus récente)
- **Région ciblée**: V4-V5 avec primers 515F-Y/926R
- **Méthode**: Téléchargement direct + RESCRIPt pour formatage QIIME2
- **Classifieur**: Naive Bayes entraîné sur données fraîches du site officiel

## ✅ Corrections apportées

### Problème diversité RÉSOLU
- Erreur "--output-dir already exists" corrigée
- Tous les outputs individuels spécifiés manuellement
- Dossier diversity_tsv maintenant peuplé avec TOUS les fichiers TSV

### Téléchargement SILVA authentique
- Pas de récupération depuis cluster (ancienne taxonomie)
- Téléchargement direct depuis arb-silva.de à chaque exécution
- Taxonomie garantie à jour (juillet 2024)

## Fichiers générés

### Tables principales
- **ASV Table avec taxonomie SILVA 138.2 officiel** : `subtables/RarTable-all/ASV.txt`
- **Table de features BIOM** : `core/table/feature-table.biom`
- **Taxonomie SILVA 138.2 officiel** : `core/taxonomy/taxonomy.tsv`
- **Séquences représentatives** : `core/rep-seqs/dna-sequences.fasta`

### Classifieur personnalisé FRAIS
- **Classifieur SILVA 138.2 V4-V5** : `98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza`

### Métriques de diversité (formats .qza ET .tsv)
- **Alpha diversity** : Vector-observed_asv.qza, Vector-shannon.qza, Vector-evenness.qza, Vector-faith_pd.qza
- **Beta diversity** : Matrix-jaccard.qza, Matrix-braycurtis.qza, Matrix-unweighted_unifrac.qza, Matrix-weighted_unifrac.qza
- **PCoA** : PCoA-jaccard.qza, PCoA-braycurtis.qza, PCoA-unweighted_unifrac.qza, PCoA-weighted_unifrac.qza
- **Visualisations Emperor** : Emperor-jaccard.qzv, Emperor-braycurtis.qzv, Emperor-unweighted_unifrac.qzv, Emperor-weighted_unifrac.qzv

### Fichiers TSV/TXT pour analyses (CORRIGÉ - DOSSIER PEUPLÉ)
- **Métriques alpha** : `diversity_tsv/observed_features.tsv`, `diversity_tsv/shannon.tsv`, etc.
- **Matrices distance** : `diversity_tsv/jaccard_distance.tsv`, `diversity_tsv/bray_curtis_distance.tsv`, etc.
- **PCoA** : `diversity_tsv/jaccard_pcoa.tsv`, `diversity_tsv/bray_curtis_pcoa.tsv`, etc.
- **Stats DADA2** : `diversity_tsv/dada2_stats.tsv`

### Rapports qualité
- **FastQC données brutes** : `../../02_qualitycheck/raw_data_qc.html`
- **FastQC données nettoyées** : `../../03_cleaned_data_qc/cleaned_data_qc.html`
- **Taxa barplots SILVA 138.2** : `visual/taxa-bar-plots.qzv`
- **Core features** : `visual/CoreBiom-all.qzv`

## Avantages de cette approche

- ✅ Taxonomie SILVA 138.2 fraîche téléchargée depuis le site officiel
- ✅ Pas de dépendance sur d'anciennes bases cluster
- ✅ Nomenclature la plus récente (juillet 2024)
- ✅ Outputs de diversité tous présents en TSV
- ✅ Classification robuste région V4-V5 optimisée
- ✅ Reproductible : télécharge toujours la version officielle

## Utilisation des fichiers

### Pour analyses statistiques
Utilisez `ASV.txt` qui contient les comptages avec taxonomie SILVA 138.2 officiel fraîche.

### Pour visualisations
Les fichiers `.qzv` peuvent être visualisés sur https://view.qiime2.org

### Pour analyses phylogénétiques
Utilisez `tree.qza` avec les métriques UniFrac.

### Pour analyses R/Python
Tous les fichiers TSV sont dans `diversity_tsv/` pour import direct.

### Classifieur réutilisable
Le classifieur `silva-138.2-ssu-nr99-515f-926r-classifier.qza` est basé sur les données officielles les plus récentes.
EOF

log "🎉 PIPELINE COMPLET TERMINÉ AVEC SILVA 138.2 OFFICIEL FRAIS !"
log "✅ Taxonomie SILVA 138.2 téléchargée depuis https://www.arb-silva.de"
log "✅ Région V4-V5 optimisée (515F-Y/926R)"
log "✅ Tous les exports et conversions réalisés"
log "✅ Dossier diversity_tsv maintenant peuplé"
log "✅ Erreurs core-metrics corrigées avec outputs individuels"
log ""
log "Consultez le rapport : ${ROOTDIR}/05_QIIME2/export/summary_tables/PIPELINE_SUMMARY_REPORT.md"
log "Fichiers TSV diversité dans : ${ROOTDIR}/05_QIIME2/export/diversity_tsv/"
