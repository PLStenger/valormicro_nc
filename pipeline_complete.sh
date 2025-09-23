#!/usr/bin/env bash

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }
log "Initialisation OK"

# ---- Fix conda issues
log "Correction des problèmes conda et Java"
if [ -f "/home/fungi/.conda/aau_token_host" ]; then
    rm -f /home/fungi/.conda/aau_token_host || log "Token non supprimé"
fi

set +u
source $(conda info --base)/etc/profile.d/conda.sh
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
set -u

# ---- 00 Génération automatique des métadonnées
log "Génération automatique des fichiers manifest et metadata"
cd "${ROOTDIR}/00_scripts"

MANIFEST="${ROOTDIR}/98_databasefiles/manifest"
MANIFEST_CONTROL="${ROOTDIR}/98_databasefiles/manifest_control"
METADATA="${ROOTDIR}/98_databasefiles/sample-metadata.tsv"

if [[ ! -f "$MANIFEST" ]] || [[ ! -f "$METADATA" ]]; then
    log "Fichiers metadata manquants - génération automatique"
    if python3 -c "import pandas" 2>/dev/null; then
        python3 "${ROOTDIR}/00_scripts/generate_qiime_files.py"
    else
        bash "${ROOTDIR}/00_scripts/generate_metadata.sh"
    fi
else
    log "Fichiers metadata existants - utilisation des fichiers présents"
fi

# ---- 01 FastQC
log "Lancement FastQC sur raw .fastq.gz"
mkdir -p "${ROOTDIR}/02_qualitycheck"
cd "${ROOTDIR}/01_raw_data"

set +u
if conda activate fastqc 2>/dev/null; then
    FASTQC_CMD="fastqc"
else
    FASTQC_CMD="conda run -n fastqc fastqc"
fi
set -u

fastq_count=0
for file in $(find . -name '*.fastq*' -type f); do
    log "FastQC: $file"
    $FASTQC_CMD "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 || continue
    fastq_count=$((fastq_count + 1))
    if [ $((fastq_count % 5)) -eq 0 ]; then sleep 2; fi
done

set +u
if conda activate multiqc 2>/dev/null; then
    cd "${ROOTDIR}/02_qualitycheck" && multiqc .
else
    cd "${ROOTDIR}/02_qualitycheck" && conda run -n multiqc multiqc .
fi
set -u

# ---- 02 Trimmomatic avec gestion correcte des fichiers paired
log "Nettoyage avec Trimmomatic - gestion des fichiers paired"
ADAPTERS="${ROOTDIR}/99_softwares/adapters/sequences.fasta"
mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

set +u
if conda activate trimmomatic 2>/dev/null; then
    TRIMMO_CMD="trimmomatic"
else
    TRIMMO_CMD="conda run -n trimmomatic trimmomatic"
fi
set -u

# Traitement avec gestion correcte des outputs Trimmomatic
pair_count=0
find . -name '*R1*.fastq*' -type f | while read r1; do
    r2="${r1/_R1/_R2}"
    if [[ -f "$r2" ]]; then
        log "Trimmomatic: pair $((++pair_count)) - $r1 et $r2"
        
        base1=$(basename "$r1")
        base2=$(basename "$r2")
        
        # Noms de sortie Trimmomatic corrects
        out1p="${ROOTDIR}/03_cleaned_data/${base1/.fastq.gz/_paired.fastq.gz}"
        out1u="${ROOTDIR}/03_cleaned_data/${base1/.fastq.gz/_unpaired.fastq.gz}"
        out2p="${ROOTDIR}/03_cleaned_data/${base2/.fastq.gz/_paired.fastq.gz}"
        out2u="${ROOTDIR}/03_cleaned_data/${base2/.fastq.gz/_unpaired.fastq.gz}"
        
        # Gestion des .fastq sans .gz
        out1p="${out1p/.fastq/_paired.fastq}"
        out1u="${out1u/.fastq/_unpaired.fastq}"
        out2p="${out2p/.fastq/_paired.fastq}"
        out2u="${out2u/.fastq/_unpaired.fastq}"
        
        $TRIMMO_CMD PE -threads 4 -phred33 "$r1" "$r2" "$out1p" "$out1u" "$out2p" "$out2u" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:150 || {
            log "ERREUR Trimmomatic sur $r1/$r2"
            continue
        }
        
        # Vérification synchronisation des fichiers paired
        log "Vérification synchronisation: $out1p et $out2p"
        if command -v seqkit >/dev/null 2>&1; then
            count1=$(seqkit stats "$out1p" | tail -n1 | awk '{print $4}')
            count2=$(seqkit stats "$out2p" | tail -n1 | awk '{print $4}')
        else
            # Méthode alternative sans seqkit
            count1=$(( $(zcat "$out1p" 2>/dev/null || cat "$out1p" | wc -l) / 4 ))
            count2=$(( $(zcat "$out2p" 2>/dev/null || cat "$out2p" | wc -l) / 4 ))
        fi
        
        log "Counts: $out1p=$count1, $out2p=$count2"
        
        if [ "$count1" != "$count2" ]; then
            log "ATTENTION: Fichiers désynchronisés détectés - tentative de correction"
            # Utiliser fastq-pair si disponible pour resynchroniser
            if command -v fastq_pair >/dev/null 2>&1; then
                log "Resynchronisation avec fastq-pair"
                fastq_pair "$out1p" "$out2p"
                # Remplacer par les fichiers synchronisés
                mv "${out1p}.paired.fq" "$out1p"
                mv "${out2p}.paired.fq" "$out2p"
                rm -f "${out1p}.single.fq" "${out2p}.single.fq"
            fi
        fi
        
        if [ $((pair_count % 3)) -eq 0 ]; then sleep 5; fi
    fi
done

# ---- Regeneration du manifest après Trimmomatic
log "IMPORTANT: Régénération du manifest avec fichiers PAIRED uniquement"
cd "${ROOTDIR}/00_scripts"

# Script pour générer un nouveau manifest avec seulement les fichiers paired
cat > regenerate_manifest_paired.py << 'EOF'
#!/usr/bin/env python3
import os
import glob
import pandas as pd

rootdir = "/nvme/bio/data_fungi/valormicro_nc"
cleaned_dir = f"{rootdir}/03_cleaned_data"
output_dir = f"{rootdir}/98_databasefiles"

# Trouver tous les fichiers paired
paired_r1 = sorted(glob.glob(f"{cleaned_dir}/*R1*paired*"))
paired_r2 = sorted(glob.glob(f"{cleaned_dir}/*R2*paired*"))

manifest_data = []
controls = []

for r1_file in paired_r1:
    # Extraire le nom d'échantillon
    basename = os.path.basename(r1_file)
    sample_id = basename.split('_')[0]  # Adapter selon votre nomenclature
    
    # Trouver le R2 correspondant
    r2_file = r1_file.replace('R1', 'R2')
    
    if os.path.exists(r2_file):
        # Vérifier si c'est un contrôle
        if any(ctrl in sample_id.lower() for ctrl in ['neg', 'blank', 'control', 'ctrl']):
            controls.append({
                'sample-id': sample_id,
                'forward-absolute-filepath': r1_file,
                'reverse-absolute-filepath': r2_file
            })
        else:
            manifest_data.append({
                'sample-id': sample_id,
                'forward-absolute-filepath': r1_file,
                'reverse-absolute-filepath': r2_file
            })

# Sauvegarder les manifests
if manifest_data:
    manifest_df = pd.DataFrame(manifest_data)
    manifest_df.to_csv(f"{output_dir}/manifest_paired", sep='\t', index=False)
    print(f"Manifest principal généré avec {len(manifest_data)} échantillons")

if controls:
    controls_df = pd.DataFrame(controls)
    controls_df.to_csv(f"{output_dir}/manifest_control_paired", sep='\t', index=False)
    print(f"Manifest contrôles généré avec {len(controls)} échantillons")
else:
    print("Aucun contrôle détecté")
EOF

python3 regenerate_manifest_paired.py

# Utiliser les nouveaux manifests
MANIFEST="${ROOTDIR}/98_databasefiles/manifest_paired"
MANIFEST_CONTROL="${ROOTDIR}/98_databasefiles/manifest_control_paired"

# ---- 03 FastQC sur données paired uniquement
log "FastQC sur données paired nettoyées"
mkdir -p "${ROOTDIR}/04_qualitycheck"
cd "${ROOTDIR}/03_cleaned_data"

set +u
if conda activate fastqc 2>/dev/null; then
    FASTQC_CMD="fastqc"
else
    FASTQC_CMD="conda run -n fastqc fastqc"
fi
set -u

# FastQC seulement sur fichiers paired
for file in $(find . -name '*paired*.fastq*' -type f); do
    log "FastQC paired: $file"
    $FASTQC_CMD "$file" -o "${ROOTDIR}/04_qualitycheck" --threads 2 || continue
done

set +u
if conda activate multiqc 2>/dev/null; then
    cd "${ROOTDIR}/04_qualitycheck" && multiqc .
else
    cd "${ROOTDIR}/04_qualitycheck" && conda run -n multiqc multiqc .
fi
set -u

# ---- 04 QIIME2 Import avec vérification
log "Importation QIIME2 avec fichiers paired synchronisés"
mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

set +u
if conda activate qiime2-2021.4 2>/dev/null; then
    QIIME_CMD=""
else
    QIIME_CMD="conda run -n qiime2-2021.4"
fi
set -u

# Vérification du manifest avant import
log "Vérification du manifest paired"
if [ ! -f "$MANIFEST" ]; then
    log "ERREUR: Manifest paired non trouvé: $MANIFEST"
    exit 1
fi

head -5 "$MANIFEST"

# Import principal avec fichiers paired
log "Import échantillons principaux (paired uniquement)"
$QIIME_CMD qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --output-path "core/demux_paired.qza" \
    --input-format PairedEndFastqManifestPhred33V2 || {
    log "ERREUR import QIIME2 paired"
    exit 1
}

# Contrôles paired si présents
HAS_CONTROLS=false
if [ -f "$MANIFEST_CONTROL" ] && [ -s "$MANIFEST_CONTROL" ]; then
    log "Import contrôles (paired uniquement)"
    $QIIME_CMD qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_CONTROL" \
        --output-path "core/demux_neg_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 && HAS_CONTROLS=true || {
        log "ERREUR import contrôles paired"
    }
fi

# ---- 05 DADA2 avec fichiers synchronisés
log "DADA2 denoising avec fichiers paired synchronisés"
cd "${ROOTDIR}/05_QIIME2/core"
$QIIME_CMD qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux_paired.qza \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads "$NTHREADS" || {
    log "ERREUR DADA2 - vérifier la synchronisation des fichiers"
    
    # Diagnostic en cas d'erreur
    log "Diagnostic des fichiers imported"
    $QIIME_CMD qiime tools export \
        --input-path demux_paired.qza \
        --output-path debug_export
    
    log "Comptage des reads dans les premiers fichiers exportés:"
    cd debug_export
    for f in $(ls *.fastq.gz | head -4); do
        count=$(( $(zcat "$f" | wc -l) / 4 ))
        echo "$f: $count reads"
    done
    
    exit 1
}

log "DADA2 réussi avec fichiers paired synchronisés!"

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
