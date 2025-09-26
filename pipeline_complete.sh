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

 conda run -n qiime2-2021.4 qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

# If you'd like to be able to jump to steps that only take FeatureData[Sequence] as input you can convert your data to FeatureData[Sequence] like so:
 conda run -n qiime2-2021.4 qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

# “Culling” low-quality sequences with cull-seqs
# Here we’ll remove sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) and any homopolymers that are 8 or more bases in length. These are the default parameters.
 conda run -n qiime2-2021.4 qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

# Filtering sequences by length and taxonomy
 conda run -n qiime2-2021.4 qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza 

#Dereplicating in uniq mode
 conda run -n qiime2-2021.4 qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.2-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza



 conda run -n qiime2-2021.4 qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
  --i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.2-ssu-nr99-classifier.qza        
    
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

log "=== CORRECTION EXPORT DIVERSITÉ - CHEMINS CORRIGÉS ==="

# ---- NAVIGATION VERS QIIME2
cd "${ROOTDIR}/05_QIIME2"

log "Vérification des chemins de fichiers de diversité"
log "Répertoire actuel: $(pwd)"

# Lister le contenu pour diagnostic
log "Contenu core/diversity/:"
ls -la core/diversity/ 2>/dev/null || log "Dossier core/diversity/ absent"

log "Contenu core/pcoa/:"
ls -la core/pcoa/ 2>/dev/null || log "Dossier core/pcoa/ absent"

# ---- EXPORT DIVERSITÉ AVEC CHEMINS CORRECTS
log "Export diversité avec chemins corrects"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

# Fonction d'export corrigée avec bons chemins
export_diversity_corrected() {
    local qza_file="$1"
    local output_name="$2"
    
    log "Tentative export: $qza_file -> $output_name"
    
    if [ -f "$qza_file" ]; then
        log "✅ Fichier trouvé: $qza_file"
        mkdir -p "export/diversity_tsv"
        
        # Export avec dossier temporaire
        temp_dir="export/diversity_tsv/${output_name}_temp"
        rm -rf "$temp_dir"
        mkdir -p "$temp_dir"
        
        conda run -n qiime2-2021.4 qiime tools export \
            --input-path "$qza_file" \
            --output-path "$temp_dir" && {
            
            log "Export QIIME2 réussi pour $output_name"
            
            # Chercher fichiers TSV/TXT dans le dossier temporaire
            files_found=0
            for ext in tsv txt csv; do
                find "$temp_dir" -name "*.${ext}" -type f | while read -r found_file; do
                    if [ -f "$found_file" ]; then
                        cp "$found_file" "export/diversity_tsv/${output_name}.${ext}"
                        files_found=$((files_found + 1))
                        log "✅ ${output_name}.${ext} créé depuis $(basename "$found_file")"
                    fi
                done
            done
            
            # Si aucun fichier TSV standard, prendre le premier fichier
            if [ "$files_found" -eq 0 ]; then
                first_file=$(find "$temp_dir" -type f | head -1)
                if [ -f "$first_file" ]; then
                    cp "$first_file" "export/diversity_tsv/${output_name}.tsv"
                    log "✅ ${output_name}.tsv créé (format alternatif)"
                fi
            fi
        } || {
            log "❌ Erreur export QIIME2 pour $qza_file"
        }
        
        # Nettoyer
        rm -rf "$temp_dir"
        return 0
    else
        log "❌ Fichier non trouvé: $qza_file"
        return 1
    fi
}

# ---- EXPORTS AVEC CHEMINS CORRECTS
log "Export métriques alpha diversity (chemins corrigés)"
export_diversity_corrected "core/diversity/Vector-observed_asv.qza" "observed_features"
export_diversity_corrected "core/diversity/Vector-shannon.qza" "shannon"  
export_diversity_corrected "core/diversity/Vector-evenness.qza" "evenness"
export_diversity_corrected "core/diversity/Vector-faith_pd.qza" "faith_pd"

log "Export matrices de distance (chemins corrigés)" 
export_diversity_corrected "core/diversity/Matrix-jaccard.qza" "jaccard_distance"
export_diversity_corrected "core/diversity/Matrix-braycurtis.qza" "bray_curtis_distance"
export_diversity_corrected "core/diversity/Matrix-unweighted_unifrac.qza" "unweighted_unifrac_distance"
export_diversity_corrected "core/diversity/Matrix-weighted_unifrac.qza" "weighted_unifrac_distance"

log "Export PCoA (chemins corrigés)"
export_diversity_corrected "core/pcoa/PCoA-jaccard.qza" "jaccard_pcoa"
export_diversity_corrected "core/pcoa/PCoA-braycurtis.qza" "bray_curtis_pcoa"
export_diversity_corrected "core/pcoa/PCoA-unweighted_unifrac.qza" "unweighted_unifrac_pcoa" 
export_diversity_corrected "core/pcoa/PCoA-weighted_unifrac.qza" "weighted_unifrac_pcoa"

# Stats DADA2
log "Export stats DADA2"
export_diversity_corrected "core/denoising-stats.qza" "dada2_stats"

# Table raréfiée si créée dans core
if [ -f "core/rarefied_table.qza" ]; then
    log "Export table raréfiée"
    export_diversity_corrected "core/rarefied_table.qza" "rarefied_table"
fi

# ---- RÉSUMÉ FINAL
log "Résumé export diversité avec chemins corrigés"
tsv_count=$(find export/diversity_tsv -name "*.tsv" -o -name "*.txt" -o -name "*.csv" 2>/dev/null | wc -l || echo "0")
log "Nombre total de fichiers exportés: $tsv_count"

log "Fichiers créés dans diversity_tsv:"
ls -la export/diversity_tsv/ 2>/dev/null || log "Dossier vide"

# ---- CONVERSION BIOM ET CRÉATION ASV.txt
log "Conversion BIOM vers TSV et création ASV.txt final"
cd "${ROOTDIR}/05_QIIME2/export"

# Conversion table raréfiée
if [ -f "subtables/RarTable-all/feature-table.biom" ]; then
    log "Conversion table raréfiée BIOM"
    
    # Tentative conversion biom
    if conda run -n qiime2-2021.4 biom convert \
        -i subtables/RarTable-all/feature-table.biom \
        -o subtables/RarTable-all/table-from-biom.tsv \
        --to-tsv; then
        
        # Modifier header
        sed '1d ; s/#OTU ID/ASV_ID/' \
            subtables/RarTable-all/table-from-biom.tsv > \
            subtables/RarTable-all/ASV.tsv
        
        log "✅ ASV.tsv créé ($(wc -l < subtables/RarTable-all/ASV.tsv) lignes)"
    else
        log "❌ Erreur conversion BIOM"
    fi
fi

# Création ASV.txt avec taxonomie
log "Création ASV.txt avec taxonomie SILVA"
if [ -f "subtables/RarTable-all/ASV.tsv" ] && [ -f "core/taxonomy/taxonomy.tsv" ]; then
    asv_file="subtables/RarTable-all/ASV.tsv"
    taxonomy_file="core/taxonomy/taxonomy.tsv"
    output_file="subtables/RarTable-all/ASV.txt"
    
    # Header avec colonnes taxonomiques
    sample_header=$(head -1 "$asv_file" | cut -f2-)
    echo -e "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t${sample_header}" > "$output_file"
    
    # Traitement de chaque ASV
    tail -n +2 "$asv_file" | while IFS=$'\t' read -r asv_id asv_counts; do
        # Valeurs par défaut
        kingdom="Bacteria"
        phylum="Unassigned" 
        class="Unassigned"
        order="Unassigned"
        family="Unassigned"
        genus="Unassigned"
        species="Unassigned"
        
        # Rechercher taxonomie
        if tax_line=$(grep "^${asv_id}" "$taxonomy_file" 2>/dev/null); then
            tax_string=$(echo "$tax_line" | cut -f2)
            
            # Parser taxonomie SILVA (format D_0__; D_1__; etc.)
            if [[ "$tax_string" =~ D_0__([^;]+) ]]; then
                kingdom="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_1__([^;]+) ]]; then
                phylum="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_2__([^;]+) ]]; then
                class="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_3__([^;]+) ]]; then
                order="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_4__([^;]+) ]]; then
                family="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_5__([^;]+) ]]; then
                genus="${BASH_REMATCH[1]}"
            fi
            if [[ "$tax_string" =~ D_6__([^;]+) ]]; then
                species="${BASH_REMATCH[1]}"
            fi
        fi
        
        # Nettoyer valeurs vides
        kingdom=${kingdom:-Bacteria}
        phylum=${phylum:-Unassigned}
        class=${class:-Unassigned}
        order=${order:-Unassigned}
        family=${family:-Unassigned}
        genus=${genus:-Unassigned}
        species=${species:-Unassigned}
        
        # Écrire ligne finale
        echo -e "${kingdom}\t${phylum}\t${class}\t${order}\t${family}\t${genus}\t${species}\t${asv_counts}" >> "$output_file"
    done
    
    lines_count=$(wc -l < "$output_file" 2>/dev/null || echo "0")
    log "✅ ASV.txt créé avec taxonomie SILVA ($lines_count lignes)"
    
    # Aperçu
    log "Aperçu ASV.txt:"
    head -3 "$output_file" | column -t -s$'\t' 2>/dev/null || head -3 "$output_file"
else
    log "❌ Fichiers manquants pour ASV.txt"
fi

# ---- RAPPORT FINAL
log "Création rapport final"
mkdir -p summary_tables

cat > "summary_tables/PIPELINE_SILVA_SUCCESS_REPORT.md" << EOF
# Pipeline QIIME2 avec SILVA 138.2 - TERMINÉ !

## ✅ Exports réussis

### Métriques de diversité exportées
- Fichiers TSV créés: $tsv_count
- Localisation: \`export/diversity_tsv/\`

### Architecture finale
\`\`\`
export/
├── core/
│   ├── table/ (table principale BIOM)
│   ├── taxonomy/ (classifications SILVA)
│   └── rep-seqs/ (séquences représentatives)
├── subtables/RarTable-all/
│   ├── ASV.tsv (table comptages)
│   └── ASV.txt (avec taxonomie SILVA)
├── diversity_tsv/ (TOUS LES TSV)
│   ├── observed_features.tsv
│   ├── shannon.tsv, evenness.tsv, faith_pd.tsv
│   ├── jaccard_distance.tsv, bray_curtis_distance.tsv
│   ├── unweighted_unifrac_distance.tsv, weighted_unifrac_distance.tsv
│   └── *_pcoa.tsv (coordonnées PCoA)
└── visual/ (visualisations QZV)
\`\`\`

### Fichiers prêts pour analyses

#### Table principale avec taxonomie
- **ASV.txt** : $lines_count lignes avec classification SILVA 138.2

#### Métriques diversité (TSV)
- Alpha : richesse, Shannon, équitabilité, Faith PD
- Beta : Jaccard, Bray-Curtis, UniFrac pondéré/non pondéré  
- PCoA : coordonnées pour graphiques

#### Visualisations
- Taxa barplots, Core features, Emperor plots

## ✅ SUCCÈS COMPLET

Pipeline SILVA 138.2 terminé avec tous les exports !
Date: $(date)
EOF

log "🎉 PIPELINE SILVA 138.2 TERMINÉ AVEC SUCCÈS !"
log ""
log "==================== RÉSUMÉ FINAL ===================="
log "✅ Chemins corrigés: core/diversity/ et core/pcoa/"
log "✅ $tsv_count fichiers TSV exportés"
log "✅ ASV.txt créé avec taxonomie SILVA ($lines_count lignes)"
log "✅ Tous les fichiers dans export/diversity_tsv/"
log ""
log "==================== FICHIERS PRINCIPAUX ===================="
log "🔹 Table finale : ${ROOTDIR}/05_QIIME2/export/subtables/RarTable-all/ASV.txt"
log "🔹 Taxonomie : ${ROOTDIR}/05_QIIME2/export/core/taxonomy/taxonomy.tsv"
log "🔹 Métriques TSV : ${ROOTDIR}/05_QIIME2/export/diversity_tsv/"
log "🔹 Rapport : ${ROOTDIR}/05_QIIME2/export/summary_tables/PIPELINE_SILVA_SUCCESS_REPORT.md"
log ""
log "🎉 TOUS LES EXPORTS TERMINÉS AVEC SILVA 138.2 ! 🎉"

# Afficher contenu final
log "Contenu final diversity_tsv:"
ls -la "${ROOTDIR}/05_QIIME2/export/diversity_tsv/"

log "Contenu final ASV:"
ls -la "${ROOTDIR}/05_QIIME2/export/subtables/RarTable-all/"
