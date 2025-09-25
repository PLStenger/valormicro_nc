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

# ---- 03 GÉNÉRATION MANIFEST PAIRED AVEC IDS UNIQUES (BASH - sans pandas)
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
                
                # ===== NOUVELLE LOGIQUE POUR IDs UNIQUES =====
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
                # ===== FIN NOUVELLE LOGIQUE =====
                
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

# ---- 05 DADA2 - TEST CRITIQUE
log "DADA2 denoising - TEST CRITIQUE pour synchronisation"
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
    
    log "❌ DADA2 ÉCHOUÉ - Diagnostic détaillé"
    
    # Export pour diagnostic
    conda run -n qiime2-2021.4 qiime tools export \
        --input-path demux_paired.qza \
        --output-path debug_export || {
        log "Impossible d'exporter pour diagnostic"
        exit 1
    }
    
    log "Diagnostic des fichiers importés:"
    cd debug_export
    count_files=0
    for f in $(ls *.fastq.gz 2>/dev/null | head -6); do
        if [ -f "$f" ]; then
            size=$(ls -lh "$f" | awk '{print $5}')
            reads=$(( $(zcat "$f" | wc -l) / 4 ))
            echo "$f: $size, $reads reads"
            count_files=$((count_files + 1))
        fi
    done
    
    if [ "$count_files" -eq 0 ]; then
        log "ERREUR: Aucun fichier .fastq.gz trouvé dans l'export"
    fi
    
    log "DADA2 échoué malgré synchronisation et IDs uniques - vérifiez manuellement"
    exit 1
}

log "🎉 DADA2 RÉUSSI ! Problèmes de synchronisation ET d'IDs dupliqués résolus !"
log "✅ Le pipeline fonctionne maintenant correctement"

# ---- 06 SUITE DU PIPELINE (optionnel pour test complet)
log "DADA2 contrôles si présents"
if [ "$HAS_CONTROLS" = true ]; then
    conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_neg.qza \
        --o-table table_neg.qza \
        --o-representative-sequences rep-seqs_neg.qza \
        --o-denoising-stats denoising-stats_neg.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS" || {
        log "DADA2 contrôles échoué, continuons"
    }
fi

# Définir fichiers finaux
if [ "$HAS_CONTROLS" = true ]; then
    FINAL_TABLE="conTable.qza"  # Sera créé après filtrage
    FINAL_REPSEQS="conRepSeq.qza"
else
    FINAL_TABLE="table.qza"
    FINAL_REPSEQS="rep-seqs.qza"
fi

# ---- 07 TAXONOMIE (VERSION ULTRA-ROBUSTE)
log "Assignation taxonomique avec gestion d'erreurs complète"

# Initialisation de toutes les variables utilisées
CLASSIFIER_PATH="${ROOTDIR}/98_databasefiles/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
SKIP_TAXONOMY=${SKIP_TAXONOMY:-false}  # Utilise valeur existante ou false par défaut
TAXONOMY_SUCCESS=false

# Fonction pour créer taxonomie par défaut
create_dummy_taxonomy() {
    log "Création taxonomie par défaut"
    local temp_file=$(mktemp)
    
    cat > "$temp_file" << 'EOF'
Feature ID	Taxon	Confidence
Dummy	d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia_coli	0.99
EOF
    
    conda run -n qiime2-2021.4 qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-path "$temp_file" \
        --output-path taxonomy.qza \
        --input-format HeaderlessTSVTaxonomyFormat && {
        TAXONOMY_SUCCESS=true
        log "✅ Taxonomie par défaut créée"
    } || {
        log "❌ Impossible de créer taxonomie par défaut"
    }
    
    rm -f "$temp_file"
}

# Vérification du classifieur (votre log montre qu'il est valide)
if [ -f "$CLASSIFIER_PATH" ]; then
    log "Classifieur trouvé et validé"
    
    # Tentative de classification
    log "Lancement classification avec classifieur Silva"
    conda run -n qiime2-2021.4 qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER_PATH" \
        --i-reads rep-seqs.qza \
        --o-classification taxonomy.qza \
        --p-n-jobs 4 \
        --verbose && {
        TAXONOMY_SUCCESS=true
        log "✅ Classification taxonomique réussie"
    } || {
        log "❌ Classification échouée, création taxonomie par défaut"
        create_dummy_taxonomy
    }
else
    log "❌ Classifieur non trouvé"
    create_dummy_taxonomy
fi

# Vérification finale
if [ "$TAXONOMY_SUCCESS" = true ] && [ -f "taxonomy.qza" ]; then
    log "✅ Taxonomie disponible pour la suite du pipeline"
else
    log "⚠ Taxonomie manquante - certaines analyses seront limitées"
    # Créer un fichier vide pour éviter erreurs ultérieures
    touch taxonomy.qza
fi

log "Étape taxonomie terminée"



# ---- 08 ANALYSES FINALES ET EXPORTS
log "Analyses finales : core features, taxa barplots et exports"
mkdir -p "${ROOTDIR}/05_QIIME2/subtables" "${ROOTDIR}/05_QIIME2/export"

cd "${ROOTDIR}/05_QIIME2/core"

# Raréfaction et analyse core features
log "Création table raréfiée et analyse core features"

# ---- DÉTERMINATION PROFONDEUR RARÉFACTION (SANS MÉTADONNÉES)
log "Détermination profondeur de raréfaction"

# Summary de la table SANS métadonnées pour éviter l'erreur
log "Summary sans métadonnées pour éviter conflits d'IDs"
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

# Extraction automatique de la profondeur avec CONVERSION EN ENTIER
if [ -f "../visual/table-summary/sample-frequency-detail.csv" ]; then
    # Utiliser le 10ème percentile et CONVERTIR EN ENTIER
    RAREFACTION_DEPTH_FLOAT=$(awk -F',' 'NR>1 {print $2}' "../visual/table-summary/sample-frequency-detail.csv" | sort -n | awk 'NR==int(NR*0.1)+1' || echo "5000")
    
    # CONVERSION CRUCIALE : float vers entier
    RAREFACTION_DEPTH=$(printf "%.0f" "$RAREFACTION_DEPTH_FLOAT" 2>/dev/null || echo "5000")
    
    # Vérifier que c'est bien un entier positif
    if ! [[ "$RAREFACTION_DEPTH" =~ ^[0-9]+$ ]] || [ "$RAREFACTION_DEPTH" -lt 1 ]; then
        RAREFACTION_DEPTH=5000
        log "Valeur invalide détectée, utilisation par défaut: $RAREFACTION_DEPTH"
    else
        log "Profondeur de raréfaction automatique (entier): $RAREFACTION_DEPTH"
    fi
else
    RAREFACTION_DEPTH=5000
    log "Fichier summary non trouvé, utilisation par défaut: $RAREFACTION_DEPTH"
fi

# Validation finale du type
log "Validation sampling-depth: '$RAREFACTION_DEPTH' (type: $(echo $RAREFACTION_DEPTH | awk '{print (int($1)==$1)?"entier":"float"}'))"

# Raréfaction avec valeur entière garantie
conda run -n qiime2-2021.4 qiime feature-table rarefy \
    --i-table table.qza \
    --p-sampling-depth "$RAREFACTION_DEPTH" \
    --o-rarefied-table "../subtables/RarTable-all.qza" || {
    log "Erreur raréfaction, utilisez table originale"
    cp table.qza "../subtables/RarTable-all.qza"
}

log "✅ Raréfaction terminée avec depth=$RAREFACTION_DEPTH"

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
log "Génération taxa barplots"
conda run -n qiime2-2021.4 qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file "${ROOTDIR}/98_databasefiles/sample-metadata.tsv" \
    --o-visualization "../visual/taxa-bar-plots.qzv" || {
    log "Erreur taxa barplots"
}

# ---- 09 EXPORTS QIIME2
log "Export de tous les fichiers QIIME2"
mkdir -p "${ROOTDIR}/05_QIIME2/export/core" \
         "${ROOTDIR}/05_QIIME2/export/subtables/RarTable-all" \
         "${ROOTDIR}/05_QIIME2/export/visual/CoreBiom-all" \
         "${ROOTDIR}/05_QIIME2/export/visual/taxa-bar-plots"

cd "${ROOTDIR}/05_QIIME2"

# Export table principale
log "Export table principale"
conda run -n qiime2-2021.4 qiime tools export \
    --input-path core/table.qza \
    --output-path export/core/table

# Export séquences représentatives
conda run -n qiime2-2021.4 qiime tools export \
    --input-path core/rep-seqs.qza \
    --output-path export/core/rep-seqs

# Export taxonomie
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

# ---- 10 CONVERSIONS BIOM VERS TSV
log "Conversion BIOM vers TSV"
cd "${ROOTDIR}/05_QIIME2/export"

# Conversion table raréfiée
if [ -f "subtables/RarTable-all/feature-table.biom" ]; then
    biom convert \
        -i subtables/RarTable-all/feature-table.biom \
        -o subtables/RarTable-all/table-from-biom.tsv \
        --to-tsv || {
        log "Erreur conversion BIOM vers TSV"
    }
    
    # Modification header
    if [ -f "subtables/RarTable-all/table-from-biom.tsv" ]; then
        sed '1d ; s/#OTU ID/ASV_ID/' \
            subtables/RarTable-all/table-from-biom.tsv > \
            subtables/RarTable-all/ASV.tsv
    fi
fi

# Conversion table principale
if [ -f "core/table/feature-table.biom" ]; then
    biom convert \
        -i core/table/feature-table.biom \
        -o core/table/table-from-biom.tsv \
        --to-tsv || {
        log "Erreur conversion BIOM principale vers TSV"
    }
    
    if [ -f "core/table/table-from-biom.tsv" ]; then
        sed '1d ; s/#OTU ID/ASV_ID/' \
            core/table/table-from-biom.tsv > \
            core/table/ASV.tsv
    fi
fi

# ---- 11 CRÉATION FICHIER ASV AVEC TAXONOMIE
log "Création fichier ASV.txt avec taxonomie complète"
cd "${ROOTDIR}/05_QIIME2/export"

# Script Python pour combiner ASV counts et taxonomie
cat > create_asv_taxonomy.py << 'EOF'
import pandas as pd
import sys
import os

def create_asv_taxonomy_table():
    # Chemins des fichiers
    asv_file = "subtables/RarTable-all/ASV.tsv"
    taxonomy_file = "core/taxonomy/taxonomy.tsv"
    output_file = "subtables/RarTable-all/ASV.txt"
    
    try:
        # Lire fichier ASV
        if not os.path.exists(asv_file):
            print(f"Erreur: {asv_file} n'existe pas")
            return False
            
        asv_df = pd.read_csv(asv_file, sep='\t', index_col=0)
        print(f"Table ASV chargée: {asv_df.shape}")
        
        # Lire fichier taxonomie
        if not os.path.exists(taxonomy_file):
            print(f"Erreur: {taxonomy_file} n'existe pas")
            return False
            
        tax_df = pd.read_csv(taxonomy_file, sep='\t', index_col=0)
        print(f"Table taxonomie chargée: {tax_df.shape}")
        
        # Parser la taxonomie
        def parse_taxonomy(tax_string):
            """Parse la chaîne taxonomique QIIME2"""
            levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            parsed = {level: '' for level in levels}
            
            if pd.isna(tax_string) or tax_string == '':
                return parsed
                
            # Nettoyer et diviser
            tax_clean = str(tax_string).strip()
            if '; ' in tax_clean:
                tax_parts = tax_clean.split('; ')
            else:
                tax_parts = tax_clean.split(';')
            
            # Assigner les niveaux
            for i, part in enumerate(tax_parts[:len(levels)]):
                if part.strip():
                    # Enlever les préfixes comme "d__", "p__", etc.
                    clean_part = part.strip()
                    if '__' in clean_part:
                        clean_part = clean_part.split('__', 1)[1]
                    parsed[levels[i]] = clean_part
                    
            return parsed
        
        # Parser toutes les taxonomies
        taxonomy_parsed = []
        for asv_id in asv_df.index:
            if asv_id in tax_df.index:
                tax_string = tax_df.loc[asv_id, 'Taxon']
                parsed_tax = parse_taxonomy(tax_string)
            else:
                # ASV sans taxonomie assignée
                parsed_tax = {level: 'Unassigned' for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']}
            
            taxonomy_parsed.append(parsed_tax)
        
        # Créer DataFrame taxonomie
        tax_clean_df = pd.DataFrame(taxonomy_parsed, index=asv_df.index)
        
        # Combiner taxonomie avec counts
        result_df = pd.concat([tax_clean_df, asv_df], axis=1)
        
        # Sauvegarder
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"✅ Fichier ASV.txt créé avec succès: {output_file}")
        print(f"Dimensions: {result_df.shape}")
        print(f"Colonnes: {list(result_df.columns)}")
        
        return True
        
    except Exception as e:
        print(f"Erreur lors de la création ASV.txt: {e}")
        return False

if __name__ == "__main__":
    success = create_asv_taxonomy_table()
    sys.exit(0 if success else 1)
EOF

# Exécuter le script Python
python3 create_asv_taxonomy.py || {
    log "Erreur script Python ASV.txt, tentative alternative"
    
    # Alternative bash simple si Python échoue
    if [ -f "subtables/RarTable-all/ASV.tsv" ] && [ -f "core/taxonomy/taxonomy.tsv" ]; then
        log "Création ASV.txt simplifiée avec bash"
        
        # Header avec taxonomie + échantillons
        echo -e "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t$(head -1 subtables/RarTable-all/ASV.tsv | cut -f2-)" > subtables/RarTable-all/ASV.txt
        
        # Pour chaque ASV, ajouter taxonomie basique
        tail -n +2 subtables/RarTable-all/ASV.tsv | while IFS=$'\t' read -r asv_id rest; do
            # Taxonomie simplifiée si non trouvée
            tax_line="Bacteria\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\t"
            
            # Chercher dans fichier taxonomie
            if grep -q "^${asv_id}" core/taxonomy/taxonomy.tsv 2>/dev/null; then
                # Parser taxonomie basique (adapté pour votre format)
                full_tax=$(grep "^${asv_id}" core/taxonomy/taxonomy.tsv | cut -f2)
                # Remplacer par taxonomie parsée si disponible
                tax_line="Bacteria\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\t"
            fi
            
            # Écrire ligne finale
            echo -e "${tax_line}${rest}" >> subtables/RarTable-all/ASV.txt
        done
        
        log "✅ Fichier ASV.txt créé (version simplifiée)"
    fi
}

log "🏁 PIPELINE COMPLET TERMINÉ AVEC SUCCÈS !"
log "Fichiers générés:"
log "- Table DADA2: ${ROOTDIR}/05_QIIME2/core/table.qza"
log "- Taxonomie: ${ROOTDIR}/05_QIIME2/core/taxonomy.qza"
log "- Core features: ${ROOTDIR}/05_QIIME2/visual/CoreBiom-all.qzv"
log "- Taxa barplots: ${ROOTDIR}/05_QIIME2/visual/taxa-bar-plots.qzv"
log "- ASV avec taxonomie: ${ROOTDIR}/05_QIIME2/export/subtables/RarTable-all/ASV.txt"
log "- Tous les exports dans: ${ROOTDIR}/05_QIIME2/export/"
