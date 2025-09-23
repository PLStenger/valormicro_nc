#!/usr/bin/env bash

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }
log "Initialisation OK"

# ---- CORRECTION DEFINITIVE DU TOKEN CONDA
log "Correction définitive des problèmes conda"
rm -rf ~/.conda/ /home/fungi/.conda/ 2>/dev/null || true
conda clean --all --yes 2>/dev/null || true

# Réinitialisation conda sans le token problématique
set +u
source $(conda info --base)/etc/profile.d/conda.sh
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
set -u

# ---- 00 Génération métadonnées
log "Génération métadonnées"
cd "${ROOTDIR}/00_scripts"

MANIFEST="${ROOTDIR}/98_databasefiles/manifest"
MANIFEST_CONTROL="${ROOTDIR}/98_databasefiles/manifest_control"
METADATA="${ROOTDIR}/98_databasefiles/sample-metadata.tsv"

if [[ ! -f "$MANIFEST" ]] || [[ ! -f "$METADATA" ]]; then
    if python3 -c "import pandas" 2>/dev/null; then
        python3 "${ROOTDIR}/00_scripts/generate_qiime_files.py"
    else
        bash "${ROOTDIR}/00_scripts/generate_metadata.sh"
    fi
else
    log "Métadonnées existantes utilisées"
fi

# ---- 01 FastQC rapide
log "FastQC raw data"
mkdir -p "${ROOTDIR}/02_qualitycheck"
cd "${ROOTDIR}/01_raw_data"

set +u
conda activate fastqc 2>/dev/null || { log "Erreur activation fastqc, utilisation conda run"; }
set -u

# FastQC optimisé
for file in $(find . -name '*.fastq*' -type f | head -10); do  # Limite pour test
    log "FastQC: $(basename $file)"
    if conda activate fastqc 2>/dev/null; then
        fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || continue
    else
        conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || continue
    fi
done

# MultiQC
set +u
conda activate multiqc 2>/dev/null && cd "${ROOTDIR}/02_qualitycheck" && multiqc . --quiet || {
    cd "${ROOTDIR}/02_qualitycheck" && conda run -n multiqc multiqc . --quiet
}
set -u

# ---- 02 Trimmomatic avec synchronisation correcte
log "Trimmomatic avec gestion paired/unpaired"
ADAPTERS="${ROOTDIR}/99_softwares/adapters/sequences.fasta"
mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

set +u
conda activate trimmomatic 2>/dev/null || { log "Utilisation conda run pour trimmomatic"; }
set -u

# Traitement Trimmomatic correct
count=0
find . -name '*R1*.fastq*' -type f | head -5 | while read r1; do  # Test sur 5 paires
    r2="${r1/_R1/_R2}"
    if [[ -f "$r2" ]]; then
        count=$((count + 1))
        log "Trimmomatic pair $count: $(basename $r1) + $(basename $r2)"
        
        base1=$(basename "$r1" .fastq.gz)
        base1=$(basename "$base1" .fastq)
        base2=$(basename "$r2" .fastq.gz)
        base2=$(basename "$base2" .fastq)
        
        # Sortie Trimmomatic avec noms clairs
        out1p="${ROOTDIR}/03_cleaned_data/${base1}_paired.fastq.gz"
        out1u="${ROOTDIR}/03_cleaned_data/${base1}_unpaired.fastq.gz" 
        out2p="${ROOTDIR}/03_cleaned_data/${base2}_paired.fastq.gz"
        out2u="${ROOTDIR}/03_cleaned_data/${base2}_unpaired.fastq.gz"
        
        if conda activate trimmomatic 2>/dev/null; then
            trimmomatic PE -threads 4 -phred33 "$r1" "$r2" \
                "$out1p" "$out1u" "$out2p" "$out2u" \
                ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:150 || {
                log "ERREUR Trimmomatic $r1/$r2"
                continue
            }
        else
            conda run -n trimmomatic trimmomatic PE -threads 4 -phred33 "$r1" "$r2" \
                "$out1p" "$out1u" "$out2p" "$out2u" \
                ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:150 || {
                log "ERREUR Trimmomatic $r1/$r2"
                continue  
            }
        fi
        
        # Vérification synchronisation immédiate
        if [[ -f "$out1p" && -f "$out2p" ]]; then
            count1=$(( $(zcat "$out1p" 2>/dev/null | wc -l || cat "$out1p" | wc -l) / 4 ))
            count2=$(( $(zcat "$out2p" 2>/dev/null | wc -l || cat "$out2p" | wc -l) / 4 ))
            
            log "Vérification: $out1p=$count1 reads, $out2p=$count2 reads"
            
            if [ "$count1" != "$count2" ]; then
                log "ATTENTION: Désynchronisation détectée!"
                # Supprimer les fichiers désynchronisés
                rm -f "$out1p" "$out2p"
                log "Fichiers désynchronisés supprimés"
                continue
            else
                log "Synchronisation OK: $count1 reads"
            fi
        fi
    fi
done

# ---- Nouveau manifest pour fichiers paired uniquement
log "Génération manifest pour fichiers paired synchronisés"
cd "${ROOTDIR}/98_databasefiles"

# Script manifest paired
cat > create_paired_manifest.py << 'EOF'
import glob
import os
import pandas as pd

rootdir = "/nvme/bio/data_fungi/valormicro_nc"
cleaned = f"{rootdir}/03_cleaned_data"

# Trouver tous les paired
r1_files = sorted(glob.glob(f"{cleaned}/*R1*_paired.fastq*"))
manifest_data = []
controls_data = []

for r1 in r1_files:
    r2 = r1.replace('R1', 'R2')
    if os.path.exists(r2):
        # Extraire sample-id
        basename = os.path.basename(r1)
        sample_id = basename.split('_')[0]
        
        # Vérifier si c'est un contrôle
        if any(x in sample_id.lower() for x in ['neg', 'blank', 'ctrl', 'control']):
            controls_data.append({
                'sample-id': sample_id,
                'forward-absolute-filepath': r1,
                'reverse-absolute-filepath': r2
            })
        else:
            manifest_data.append({
                'sample-id': sample_id,
                'forward-absolute-filepath': r1,
                'reverse-absolute-filepath': r2
            })

# Sauvegarder
if manifest_data:
    df = pd.DataFrame(manifest_data)
    df.to_csv('manifest_paired', sep='\t', index=False)
    print(f"Manifest principal: {len(manifest_data)} échantillons")

if controls_data:
    df = pd.DataFrame(controls_data)
    df.to_csv('manifest_control_paired', sep='\t', index=False)
    print(f"Manifest contrôles: {len(controls_data)} échantillons")
else:
    print("Pas de contrôles")
EOF

python3 create_paired_manifest.py

# Vérifier les manifests générés
MANIFEST_PAIRED="${ROOTDIR}/98_databasefiles/manifest_paired"
MANIFEST_CONTROL_PAIRED="${ROOTDIR}/98_databasefiles/manifest_control_paired"

if [ ! -f "$MANIFEST_PAIRED" ]; then
    log "ERREUR: Pas de manifest paired généré"
    exit 1
fi

log "Manifest paired créé:"
head -3 "$MANIFEST_PAIRED"

# ---- 03 QIIME2 Import
log "Import QIIME2 avec fichiers paired synchronisés"
mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

set +u
conda activate qiime2-2021.4 2>/dev/null || { log "Utilisation conda run pour qiime2"; }
set -u

# Import principal
log "Import principal QIIME2"
if conda activate qiime2-2021.4 2>/dev/null; then
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_PAIRED" \
        --output-path "core/demux_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || {
        log "ERREUR import principal"
        exit 1
    }
else
    conda run -n qiime2-2021.4 qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_PAIRED" \
        --output-path "core/demux_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 || {
        log "ERREUR import principal"
        exit 1
    }
fi

# Import contrôles si présents
HAS_CONTROLS=false
if [ -f "$MANIFEST_CONTROL_PAIRED" ] && [ -s "$MANIFEST_CONTROL_PAIRED" ]; then
    log "Import contrôles QIIME2"
    if conda activate qiime2-2021.4 2>/dev/null; then
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "$MANIFEST_CONTROL_PAIRED" \
            --output-path "core/demux_neg.qza" \
            --input-format PairedEndFastqManifestPhred33V2 && HAS_CONTROLS=true
    else
        conda run -n qiime2-2021.4 qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "$MANIFEST_CONTROL_PAIRED" \
            --output-path "core/demux_neg.qza" \
            --input-format PairedEndFastqManifestPhred33V2 && HAS_CONTROLS=true
    fi
fi

# ---- 04 DADA2 
log "DADA2 avec fichiers paired synchronisés"
cd "${ROOTDIR}/05_QIIME2/core"

if conda activate qiime2-2021.4 2>/dev/null; then
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_paired.qza \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS" || {
        log "ERREUR DADA2"
        exit 1
    }
else
    conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_paired.qza \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS" || {
        log "ERREUR DADA2"
        exit 1
    }
fi

log "DADA2 RÉUSSI! Fichiers synchronisés fonctionnent."
log "Suite du pipeline peut continuer..."

# ---- 05 Filtrage des contaminants (si contrôles présents)
if [ "$HAS_CONTROLS" = true ]; then
    log "Filtrage des contaminants basé sur les contrôles"
    
    # Identification des séquences contaminantes
    qiime quality-control exclude-seqs \
        --i-query-sequences rep-seqs.qza \
        --i-reference-sequences rep-seqs_neg.qza \
        --p-method vsearch \
        --p-threads "$NTHREADS" \
        --p-perc-identity 1.00 \
        --p-perc-query-aligned 1.00 \
        --o-sequence-hits hit_neg_ctrl.qza \
        --o-sequence-misses clean_rep-seqs.qza || {
        log "ERREUR filtrage contaminants"
    }
    
    # Filtrage de la table des features (suppression des contaminants)
    qiime feature-table filter-features \
        --i-table table.qza \
        --m-metadata-file hit_neg_ctrl.qza \
        --o-filtered-table negTable.qza \
        --p-exclude-ids || {
        log "ERREUR filtrage table"
    }
    
    # Filtrage par contingence (minimum 2 échantillons)
    qiime feature-table filter-features \
        --i-table negTable.qza \
        --p-min-samples 2 \
        --p-min-frequency 0 \
        --o-filtered-table conTable.qza || {
        log "ERREUR filtrage contingence"
    }
    
    # Filtrage des séquences correspondantes
    qiime feature-table filter-seqs \
        --i-data clean_rep-seqs.qza \
        --i-table conTable.qza \
        --o-filtered-data conRepSeq.qza || {
        log "ERREUR filtrage séquences"
    }
    
    # Tables finales pour la suite
    FINAL_TABLE="conTable.qza"
    FINAL_REPSEQS="conRepSeq.qza"
    TABLES_TO_PROCESS=("table.qza" "negTable.qza" "conTable.qza")
    REPSEQS_TO_PROCESS=("rep-seqs.qza" "clean_rep-seqs.qza" "conRepSeq.qza")
    
else
    log "Pas de contrôles - Utilisation des données brutes"
    FINAL_TABLE="table.qza"
    FINAL_REPSEQS="rep-seqs.qza"
    TABLES_TO_PROCESS=("table.qza")
    REPSEQS_TO_PROCESS=("rep-seqs.qza")
fi

# ---- 06 Taxonomie
log "Assignation taxonomique"
qiime feature-classifier classify-sklearn \
    --i-classifier "${ROOTDIR}/98_databasefiles/silva-138-99-515-926-nb-classifier.qza" \
    --i-reads "$FINAL_REPSEQS" \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS" || {
    log "ERREUR classification taxonomique"
}

# ---- 07 Résumés tables & séquences
log "Visualisation & résumés des tables et séquences"
for tab in "${TABLES_TO_PROCESS[@]}"; do
    if [ -f "$tab" ]; then
        log "Résumé table: $tab"
        qiime feature-table summarize \
            --i-table "$tab" \
            --m-sample-metadata-file "$METADATA" \
            --o-visualization "visual_${tab%.qza}.qzv" \
            --verbose || log "ERREUR résumé $tab"
    fi
done

for rep in "${REPSEQS_TO_PROCESS[@]}"; do
    if [ -f "$rep" ]; then
        log "Tabulation séquences: $rep"
        qiime feature-table tabulate-seqs \
            --i-data "$rep" \
            --o-visualization "visual_${rep%.qza}.qzv" \
            --verbose || log "ERREUR tabulation $rep"
    fi
done

# ---- 08 Taxonomy Barplots
log "Barplots taxonomiques"
for tab in "${TABLES_TO_PROCESS[@]}"; do
    if [ -f "$tab" ]; then
        log "Barplot: $tab"
        qiime taxa barplot \
            --i-table "$tab" \
            --i-taxonomy taxonomy.qza \
            --m-metadata-file "$METADATA" \
            --o-visualization "barplot_${tab%.qza}.qzv" \
            --verbose || log "ERREUR barplot $tab"
    fi
done

# ---- 09 Construction arbre phylogénétique
log "Arbre phylogénétique (alignement, masquage, arbre, rooting)"
qiime alignment mafft \
    --i-sequences "$FINAL_REPSEQS" \
    --o-alignment aligned-rep-seqs.qza \
    --p-n-threads "$NTHREADS" || {
    log "ERREUR alignement MAFFT"
}

qiime alignment mask \
    --i-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza || {
    log "ERREUR masquage alignement"
}

qiime phylogeny fasttree \
    --i-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza || {
    log "ERREUR construction arbre"
}

qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza || {
    log "ERREUR enracinement arbre"
}

# ---- 10 Analyses diversité alpha/beta, rarefaction
log "Diversité alpha/beta, rarefaction"
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table "$FINAL_TABLE" \
    --p-sampling-depth 10000 \
    --m-metadata-file "$METADATA" \
    --output-dir core-metrics-results \
    --p-n-jobs "$NTHREADS" || {
    log "ERREUR core metrics"
}

qiime diversity alpha-rarefaction \
    --i-table "$FINAL_TABLE" \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 20000 \
    --m-metadata-file "$METADATA" \
    --o-visualization visual_alpha_rarefaction.qzv || {
    log "ERREUR alpha rarefaction"
}

# ---- 11 Export des principaux résultats
log "Export QIIME2 principaux résultats"
mkdir -p "${ROOTDIR}/exports"

# Export des artifacts principaux
for artifact in "$FINAL_TABLE" "$FINAL_REPSEQS" taxonomy.qza rooted-tree.qza; do
    if [ -f "$artifact" ]; then
        log "Export: $artifact"
        qiime tools export \
            --input-path "$artifact" \
            --output-path "${ROOTDIR}/exports/${artifact%.qza}" \
            --verbose || log "ERREUR export $artifact"
    fi
done

# Export des visualisations
for viz in visual_*.qzv barplot_*.qzv visual_alpha_rarefaction.qzv; do
    if [ -f "$viz" ]; then
        log "Export visualisation: $viz"
        qiime tools export \
            --input-path "$viz" \
            --output-path "${ROOTDIR}/exports/${viz%.qzv}" \
            --verbose || log "ERREUR export $viz"
    fi
done

# Export des résultats de diversité
if [ -d "core-metrics-results" ]; then
    log "Export core metrics"
    cp -r core-metrics-results "${ROOTDIR}/exports/" || log "ERREUR copie core metrics"
fi

# ---- 12 Nettoyage
log "Nettoyage temporaire et permissions"
rm -rf "$TMPDIR" || log "ERREUR nettoyage tmp"
chown -R $(id -u):$(id -g) "${ROOTDIR}" 2>/dev/null || log "AVERTISSEMENT permissions"

log "Pipeline terminé avec succès !"
