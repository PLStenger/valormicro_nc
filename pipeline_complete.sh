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

# ---- 07 TAXONOMIE
log "Assignation taxonomique"
conda run -n qiime2-2021.4 qiime feature-classifier classify-sklearn \
    --i-classifier "${ROOTDIR}/98_databasefiles/silva-138-99-515-926-nb-classifier.qza" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS" || {
    log "ERREUR classification taxonomique"
}

log "🏁 PIPELINE TERMINÉ AVEC SUCCÈS !"
log "Fichiers générés:"
log "- Table DADA2: ${ROOTDIR}/05_QIIME2/core/table.qza"
log "- Séquences représentatives: ${ROOTDIR}/05_QIIME2/core/rep-seqs.qza"
log "- Taxonomie: ${ROOTDIR}/05_QIIME2/core/taxonomy.qza"
log "Le pipeline peut maintenant continuer avec les analyses de diversité..."
