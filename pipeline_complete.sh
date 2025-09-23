#!/usr/bin/env bash
#SBATCH --job-name=valormicro
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=local
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/nvme/bio/data_fungi/valormicro_nc/00_scripts/valormicro.err"
#SBATCH --output="/nvme/bio/data_fungi/valormicro_nc/00_scripts/valormicro.out"

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
export NTHREADS=36   # Utilise tous les CPUs alloués
export TMPDIR="${ROOTDIR}/tmp"
mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }

# ---- 00 Génération automatique des métadonnées
log "Génération automatique des fichiers manifest et metadata"
cd "${ROOTDIR}/00_scripts"

# Vérifier si les fichiers existent déjà
MANIFEST="${ROOTDIR}/98_databasefiles/manifest"
MANIFEST_CONTROL="${ROOTDIR}/98_databasefiles/manifest_control"
METADATA="${ROOTDIR}/98_databasefiles/sample-metadata.tsv"

if [[ ! -f "$MANIFEST" ]] || [[ ! -f "$METADATA" ]]; then
    log "Fichiers metadata manquants - génération automatique"
    
    # Utiliser le script Python si pandas est disponible, sinon bash
    if python3 -c "import pandas" 2>/dev/null; then
        log "Utilisation du générateur Python"
        python3 "${ROOTDIR}/00_scripts/generate_qiime_files.py"
    else
        log "Utilisation du générateur Bash (pandas non disponible)"
        bash "${ROOTDIR}/00_scripts/generate_metadata.sh"
    fi
else
    log "Fichiers metadata existants - utilisation des fichiers présents"
fi

# ---- 01 FastQC
log "Lancement FastQC sur raw .fastq.gz"
mkdir -p "${ROOTDIR}/02_qualitycheck"
cd "${ROOTDIR}/01_raw_data"
find . -name '*.fastq*' -type f | parallel -j "$NTHREADS" "conda run -n fastqc fastqc {} -o ${ROOTDIR}/02_qualitycheck"

log "MultiQC - resume FastQC"
cd "${ROOTDIR}/02_qualitycheck"
conda run -n multiqc multiqc .

# ---- 02 Trimmomatic
log "Nettoyage avec Trimmomatic"
ADAPTERS="${ROOTDIR}/99_softwares/adapters/sequences.fasta"
mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

# Fonction pour traiter un pair de fichiers
process_pair() {
    local r1="$1"
    local r2="$2"
    local out1p="${ROOTDIR}/03_cleaned_data/${r1##*/}"
    local out1u="${ROOTDIR}/03_cleaned_data/${r1##*/}"
    local out2p="${ROOTDIR}/03_cleaned_data/${r2##*/}"
    local out2u="${ROOTDIR}/03_cleaned_data/${r2##*/}"
    
    out1p="${out1p/.fastq/_paired.fastq}"
    out1u="${out1u/.fastq/_unpaired.fastq}"
    out2p="${out2p/.fastq/_paired.fastq}"
    out2u="${out2u/.fastq/_unpaired.fastq}"
    
    conda run -n trimmomatic trimmomatic PE -threads 4 -phred33 "$r1" "$r2" "$out1p" "$out1u" "$out2p" "$out2u" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:2:30 MINLEN:150
}

export -f process_pair
export ROOTDIR ADAPTERS

# Traitement parallélisé des pairs
find . -name '*R1*.fastq*' -type f | while read r1; do
    r2="${r1/_R1/_R2}"
    if [[ -f "$r2" ]]; then
        echo "$r1 $r2"
    fi
done | parallel -j $(( NTHREADS / 4 )) --colsep ' ' process_pair {1} {2}

log "FastQC + MultiQC sur reads nettoyés"
mkdir -p "${ROOTDIR}/04_qualitycheck"
cd "${ROOTDIR}/03_cleaned_data"
find . -name '*.fastq*' -type f | parallel -j "$NTHREADS" "conda run -n fastqc fastqc {} -o ${ROOTDIR}/04_qualitycheck"
cd "${ROOTDIR}/04_qualitycheck"
conda run -n multiqc multiqc .

# ---- 03 QIIME2 Import - Détection automatique des contrôles
log "Importation dans QIIME2 avec détection automatique des contrôles"
mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

# Import des échantillons principaux
conda run -n qiime2-2021.4 qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --output-path "core/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2

# Vérification et import des contrôles si présents
HAS_CONTROLS=false
if [ -f "$MANIFEST_CONTROL" ] && [ -s "$MANIFEST_CONTROL" ]; then
    log "Contrôles détectés - Import des échantillons de contrôle"
    conda run -n qiime2-2021.4 qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_CONTROL" \
        --output-path "core/demux_neg.qza" \
        --input-format PairedEndFastqManifestPhred33V2
    HAS_CONTROLS=true
else
    log "Aucun contrôle détecté - Poursuite sans contrôles"
fi

# ---- 04 QIIME2 Denoising (DADA2)
log "DADA2 denoise - échantillons principaux"
cd "${ROOTDIR}/05_QIIME2/core"
conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "$NTHREADS"

# Denoising des contrôles si présents
if [ "$HAS_CONTROLS" = true ]; then
    log "DADA2 denoise - échantillons de contrôle"
    conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_neg.qza \
        --o-table table_neg.qza \
        --o-representative-sequences rep-seqs_neg.qza \
        --o-denoising-stats denoising-stats_neg.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS"
fi

# ---- 05 Filtrage des contaminants (si contrôles présents)
if [ "$HAS_CONTROLS" = true ]; then
    log "Filtrage des contaminants basé sur les contrôles"
    
    # Identification des séquences contaminantes
    conda run -n qiime2-2021.4 qiime quality-control exclude-seqs \
        --i-query-sequences rep-seqs.qza \
        --i-reference-sequences rep-seqs_neg.qza \
        --p-method vsearch \
        --p-threads "$NTHREADS" \
        --p-perc-identity 1.00 \
        --p-perc-query-aligned 1.00 \
        --o-sequence-hits hit_neg_ctrl.qza \
        --o-sequence-misses clean_rep-seqs.qza
    
    # Filtrage de la table des features (suppression des contaminants)
    conda run -n qiime2-2021.4 qiime feature-table filter-features \
        --i-table table.qza \
        --m-metadata-file hit_neg_ctrl.qza \
        --o-filtered-table negTable.qza \
        --p-exclude-ids
    
    # Filtrage par contingence (minimum 2 échantillons)
    conda run -n qiime2-2021.4 qiime feature-table filter-features \
        --i-table negTable.qza \
        --p-min-samples 2 \
        --p-min-frequency 0 \
        --o-filtered-table conTable.qza
    
    # Filtrage des séquences correspondantes
    conda run -n qiime2-2021.4 qiime feature-table filter-seqs \
        --i-data clean_rep-seqs.qza \
        --i-table conTable.qza \
        --o-filtered-data conRepSeq.qza
    
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
conda run -n qiime2-2021.4 qiime feature-classifier classify-sklearn \
    --i-classifier "${ROOTDIR}/98_databasefiles/silva-138-99-515-926-nb-classifier.qza" \
    --i-reads "$FINAL_REPSEQS" \
    --o-classification taxonomy.qza \
    --p-n-jobs "$NTHREADS"

# ---- 07 Résumés tables & séquences
log "Visualisation & résumés des tables et séquences"
for tab in "${TABLES_TO_PROCESS[@]}"; do
    if [ -f "$tab" ]; then
        conda run -n qiime2-2021.4 qiime feature-table summarize \
            --i-table "$tab" \
            --m-sample-metadata-file "$METADATA" \
            --o-visualization "visual_${tab%.qza}.qzv" \
            --verbose
    fi
done

for rep in "${REPSEQS_TO_PROCESS[@]}"; do
    if [ -f "$rep" ]; then
        conda run -n qiime2-2021.4 qiime feature-table tabulate-seqs \
            --i-data "$rep" \
            --o-visualization "visual_${rep%.qza}.qzv" \
            --verbose
    fi
done

# ---- 08 Taxonomy Barplots
log "Barplots taxonomiques"
for tab in "${TABLES_TO_PROCESS[@]}"; do
    if [ -f "$tab" ]; then
        conda run -n qiime2-2021.4 qiime taxa barplot \
            --i-table "$tab" \
            --i-taxonomy taxonomy.qza \
            --m-metadata-file "$METADATA" \
            --o-visualization "barplot_${tab%.qza}.qzv" \
            --verbose
    fi
done

# ---- 09 Construction arbre phylogénétique
log "Arbre phylogénétique (alignement, masquage, arbre, rooting)"
conda run -n qiime2-2021.4 qiime alignment mafft \
    --i-sequences "$FINAL_REPSEQS" \
    --o-alignment aligned-rep-seqs.qza \
    --p-n-threads "$NTHREADS"

conda run -n qiime2-2021.4 qiime alignment mask \
    --i-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza

conda run -n qiime2-2021.4 qiime phylogeny fasttree \
    --i-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza

conda run -n qiime2-2021.4 qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

# ---- 10 Analyses diversité alpha/beta, rarefaction
log "Diversité alpha/beta, rarefaction"
conda run -n qiime2-2021.4 qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table "$FINAL_TABLE" \
    --p-sampling-depth 10000 \
    --m-metadata-file "$METADATA" \
    --output-dir core-metrics-results \
    --p-n-jobs "$NTHREADS"

conda run -n qiime2-2021.4 qiime diversity alpha-rarefaction \
    --i-table "$FINAL_TABLE" \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 20000 \
    --m-metadata-file "$METADATA" \
    --o-visualization visual_alpha_rarefaction.qzv

# ---- 11 Export des principaux résultats
log "Export QIIME2 principaux résultats"
mkdir -p "${ROOTDIR}/exports"

# Export des artifacts principaux
for artifact in "$FINAL_TABLE" "$FINAL_REPSEQS" taxonomy.qza rooted-tree.qza; do
    if [ -f "$artifact" ]; then
        conda run -n qiime2-2021.4 qiime tools export \
            --input-path "$artifact" \
            --output-path "${ROOTDIR}/exports/${artifact%.qza}" \
            --verbose
    fi
done

# Export des visualisations
for viz in visual_*.qzv barplot_*.qzv visual_alpha_rarefaction.qzv; do
    if [ -f "$viz" ]; then
        conda run -n qiime2-2021.4 qiime tools export \
            --input-path "$viz" \
            --output-path "${ROOTDIR}/exports/${viz%.qzv}" \
            --verbose
    fi
done

# Export des résultats de diversité
if [ -d "core-metrics-results" ]; then
    cp -r core-metrics-results "${ROOTDIR}/exports/"
fi

# ---- 12 Nettoyage
log "Nettoyage temporaire et permissions"
rm -rf "$TMPDIR"
chown -R $(id -u):$(id -g) "${ROOTDIR}" 2>/dev/null || true

log "Pipeline terminé avec succès !"
