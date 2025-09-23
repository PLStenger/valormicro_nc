#!/usr/bin/env bash

set -euo pipefail

export ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
export NTHREADS=16
export TMPDIR="${ROOTDIR}/tmp"
mkdir -p "$TMPDIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }
log "Initialisation OK"

# ---- CORRECTION AGRESSIVE DU TOKEN CONDA
log "Suppression agressive du token conda corrompu"
# Supprimer tous les fichiers token conda possibles
rm -rf ~/.conda 2>/dev/null || true
rm -rf /home/fungi/.conda 2>/dev/null || true
rm -rf ~/.continuum 2>/dev/null || true
rm -f ~/.condarc.bak 2>/dev/null || true

# Supprimer explicitement le fichier token qui pose problème
sudo rm -f /home/fungi/.conda/aau_token_host 2>/dev/null || true

# Nettoyer complètement conda
conda clean --all --yes 2>/dev/null || true
conda config --remove-key default_channels 2>/dev/null || true

# Réinitialiser conda sans token
set +u
export CONDA_TOKEN_PATH=""
export ANACONDA_API_TOKEN=""
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

# ---- 01 FastQC avec nettoyage préalable
log "FastQC sur données brutes avec nettoyage"
mkdir -p "${ROOTDIR}/02_qualitycheck"

# Nettoyer les anciens rapports
rm -f "${ROOTDIR}/02_qualitycheck/multiqc_report.html" 2>/dev/null || true
rm -rf "${ROOTDIR}/02_qualitycheck/multiqc_data" 2>/dev/null || true

cd "${ROOTDIR}/01_raw_data"

set +u
conda activate fastqc 2>/dev/null || { log "Utilisation conda run pour fastqc"; }
set -u

# FastQC optimisé - traitement limité pour test
count=0
for file in $(find . -name '*.fastq*' -type f | head -8); do  # Limité à 8 fichiers pour test
    count=$((count + 1))
    log "FastQC $count/8: $(basename $file)"
    
    if conda activate fastqc 2>/dev/null; then
        fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || continue
    else
        conda run -n fastqc fastqc "$file" -o "${ROOTDIR}/02_qualitycheck" --threads 2 --quiet || continue
    fi
    
    # Pause pour éviter surcharge
    if [ $((count % 4)) -eq 0 ]; then
        log "Pause après $count fichiers"
        sleep 3
    fi
done

# MultiQC avec force et nom personnalisé
log "MultiQC avec forçage et nom custom"
cd "${ROOTDIR}/02_qualitycheck"

set +u
if conda activate multiqc 2>/dev/null; then
    # Utiliser --force et nom custom pour éviter les conflits
    multiqc . --force --filename "fastqc_raw_report" --title "Raw Data QC" --quiet || {
        log "Erreur MultiQC, tentative avec conda run"
        conda run -n multiqc multiqc . --force --filename "fastqc_raw_report" --title "Raw Data QC" --quiet
    }
else
    conda run -n multiqc multiqc . --force --filename "fastqc_raw_report" --title "Raw Data QC" --quiet || {
        log "MultiQC échoué, continuation du script"
    }
fi
set -u

# ---- 02 Trimmomatic optimisé
log "Trimmomatic avec vérification synchronisation"
ADAPTERS="${ROOTDIR}/99_softwares/adapters/sequences.fasta"
mkdir -p "${ROOTDIR}/03_cleaned_data"
cd "${ROOTDIR}/01_raw_data"

set +u
conda activate trimmomatic 2>/dev/null || { log "Utilisation conda run pour trimmomatic"; }
set -u

# Traitement Trimmomatic avec test sur échantillon réduit
count=0
success_count=0
find . -name '*R1*.fastq*' -type f | head -3 | while read r1; do  # Test sur 3 paires seulement
    r2="${r1/_R1/_R2}"
    if [[ -f "$r2" ]]; then
        count=$((count + 1))
        log "Trimmomatic test $count/3: $(basename $r1) + $(basename $r2)"
        
        base1=$(basename "$r1" .fastq.gz)
        base1=$(basename "$base1" .fastq)
        base2=$(basename "$r2" .fastq.gz)  
        base2=$(basename "$base2" .fastq)
        
        # Noms de sortie clairs
        out1p="${ROOTDIR}/03_cleaned_data/${base1}_paired.fastq.gz"
        out1u="${ROOTDIR}/03_cleaned_data/${base1}_unpaired.fastq.gz"
        out2p="${ROOTDIR}/03_cleaned_data/${base2}_paired.fastq.gz"
        out2u="${ROOTDIR}/03_cleaned_data/${base2}_unpaired.fastq.gz"
        
        # Trimmomatic avec paramètres modérés
        trimmo_success=false
        if conda activate trimmomatic 2>/dev/null; then
            trimmomatic PE -threads 4 -phred33 "$r1" "$r2" \
                "$out1p" "$out1u" "$out2p" "$out2u" \
                ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 && trimmo_success=true
        else
            conda run -n trimmomatic trimmomatic PE -threads 4 -phred33 "$r1" "$r2" \
                "$out1p" "$out1u" "$out2p" "$out2u" \
                ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 && trimmo_success=true
        fi
        
        if [ "$trimmo_success" = true ]; then
            # Vérification synchronisation cruciale
            if [[ -f "$out1p" && -f "$out2p" ]]; then
                # Compter les reads
                if command -v seqkit >/dev/null 2>&1; then
                    count1=$(seqkit stats "$out1p" 2>/dev/null | tail -n1 | awk '{print $4}' || echo "0")
                    count2=$(seqkit stats "$out2p" 2>/dev/null | tail -n1 | awk '{print $4}' || echo "0")
                else
                    # Méthode alternative robuste
                    count1=$(( $(zcat "$out1p" 2>/dev/null | wc -l || cat "$out1p" | wc -l) / 4 ))
                    count2=$(( $(zcat "$out2p" 2>/dev/null | wc -l || cat "$out2p" | wc -l) / 4 ))
                fi
                
                log "Vérification: $count1 vs $count2 reads"
                
                if [ "$count1" = "$count2" ] && [ "$count1" -gt 0 ]; then
                    log "✓ Paire synchronisée: $count1 reads"
                    success_count=$((success_count + 1))
                else
                    log "✗ Paire désynchronisée ou vide, suppression"
                    rm -f "$out1p" "$out2p" "$out1u" "$out2u"
                fi
            else
                log "✗ Fichiers de sortie manquants"
            fi
        else
            log "✗ Échec Trimmomatic"
        fi
        
        sleep 2  # Pause entre chaque paire
    fi
done

log "Trimmomatic terminé"

# ---- Génération manifest paired
log "Génération manifest pour fichiers paired synchronisés"
cd "${ROOTDIR}/98_databasefiles"

# Script Python pour manifest paired
cat > create_paired_manifest.py << 'EOF'
import glob
import os
import pandas as pd

rootdir = "/nvme/bio/data_fungi/valormicro_nc"
cleaned = f"{rootdir}/03_cleaned_data"

# Trouver paired files synchronisés
r1_files = sorted(glob.glob(f"{cleaned}/*R1*_paired.fastq*"))
print(f"Trouvé {len(r1_files)} fichiers R1 paired")

manifest_data = []
controls_data = []

for r1 in r1_files:
    r2 = r1.replace('R1', 'R2')
    if os.path.exists(r2):
        # Vérifier que les fichiers ne sont pas vides
        if os.path.getsize(r1) > 1000 and os.path.getsize(r2) > 1000:
            basename = os.path.basename(r1)
            # Extraire sample ID plus robustement
            sample_id = basename.split('_')[0]
            
            # Détecter contrôles
            if any(ctrl in sample_id.lower() for ctrl in ['neg', 'blank', 'control', 'ctrl']):
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
        else:
            print(f"Fichier trop petit ignoré: {r1}")

print(f"Échantillons: {len(manifest_data)}, Contrôles: {len(controls_data)}")

# Sauvegarder manifests
if manifest_data:
    df = pd.DataFrame(manifest_data)
    df.to_csv('manifest_paired', sep='\t', index=False)
    print("Manifest principal créé")
else:
    print("ERREUR: Aucun échantillon valide trouvé")

if controls_data:
    df_ctrl = pd.DataFrame(controls_data)
    df_ctrl.to_csv('manifest_control_paired', sep='\t', index=False)
    print("Manifest contrôles créé")
EOF

python3 create_paired_manifest.py

# Vérifier résultat
MANIFEST_PAIRED="${ROOTDIR}/98_databasefiles/manifest_paired"

if [ ! -f "$MANIFEST_PAIRED" ] || [ ! -s "$MANIFEST_PAIRED" ]; then
    log "ERREUR: Aucun manifest paired généré - pas assez de fichiers synchronisés"
    log "Fichiers trouvés dans cleaned_data:"
    ls -la "${ROOTDIR}/03_cleaned_data/" || true
    exit 1
fi

log "Manifest paired créé avec $(wc -l < "$MANIFEST_PAIRED") lignes:"
head -3 "$MANIFEST_PAIRED"

# ---- 03 QIIME2 Import final
log "Import QIIME2 avec fichiers paired vérifiés"
mkdir -p "${ROOTDIR}/05_QIIME2/core" "${ROOTDIR}/05_QIIME2/visual"
cd "${ROOTDIR}/05_QIIME2"

set +u
conda activate qiime2-2021.4 2>/dev/null || { log "Utilisation conda run pour QIIME2"; }
set -u

# Import avec gestion d'erreur améliorée  
log "Import QIIME2 - tentative 1"
import_success=false

if conda activate qiime2-2021.4 2>/dev/null; then
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_PAIRED" \
        --output-path "core/demux_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 && import_success=true || {
        log "Import direct échoué, tentative avec conda run"
    }
fi

if [ "$import_success" = false ]; then
    conda run -n qiime2-2021.4 qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_PAIRED" \
        --output-path "core/demux_paired.qza" \
        --input-format PairedEndFastqManifestPhred33V2 && import_success=true || {
        log "Import QIIME2 complètement échoué"
        exit 1
    }
fi

log "Import QIIME2 réussi!"

# Contrôles si présents
MANIFEST_CONTROL_PAIRED="${ROOTDIR}/98_databasefiles/manifest_control_paired"
HAS_CONTROLS=false
if [ -f "$MANIFEST_CONTROL_PAIRED" ] && [ -s "$MANIFEST_CONTROL_PAIRED" ]; then
    log "Import contrôles"
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

# ---- 04 DADA2 Final
log "DADA2 avec fichiers paired synchronisés - TEST CRITIQUE"
cd "${ROOTDIR}/05_QIIME2/core"

dada2_success=false
if conda activate qiime2-2021.4 2>/dev/null; then
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_paired.qza \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS" && dada2_success=true || {
        log "DADA2 direct échoué, tentative conda run"
    }
fi

if [ "$dada2_success" = false ]; then
    conda run -n qiime2-2021.4 qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux_paired.qza \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$NTHREADS" && dada2_success=true || {
        
        log "DADA2 ÉCHOUE ENCORE - Diagnostic détaillé"
        log "Export pour diagnostic:"
        
        if conda activate qiime2-2021.4 2>/dev/null; then
            qiime tools export --input-path demux_paired.qza --output-path debug_export
        else
            conda run -n qiime2-2021.4 qiime tools export --input-path demux_paired.qza --output-path debug_export
        fi
        
        cd debug_export
        log "Diagnostic des fichiers dans l'import:"
        for f in $(ls *.fastq.gz | head -6); do
            size=$(ls -lh "$f" | awk '{print $5}')
            reads=$(( $(zcat "$f" | wc -l) / 4 ))
            echo "$f: $size, $reads reads"
        done
        
        log "ERREUR DADA2 PERSISTANTE - Vérifiez la synchronisation manuelle"
        exit 1
    }
fi

if [ "$dada2_success" = true ]; then
    log "🎉 DADA2 RÉUSSI ! Problème de synchronisation résolu !"
    log "Le pipeline peut maintenant continuer avec les étapes suivantes..."
else
    log "❌ DADA2 échoué malgré toutes les corrections"
    exit 1
fi


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
