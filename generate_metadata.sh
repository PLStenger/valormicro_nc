#!/usr/bin/env bash

set -euo pipefail

# Configuration
ROOTDIR="/nvme/bio/data_fungi/valormicro_nc"
RAW_DATA_DIR="${ROOTDIR}/01_raw_data"
DATABASE_DIR="${ROOTDIR}/98_databasefiles"
mkdir -p "$DATABASE_DIR"

log() { echo -e "\n[$(date +'%F %T')] $*\n"; }

# Function pour extraire le nom d'échantillon depuis le nom de fichier
extract_sample_name() {
    local filename="$1"
    # Pour les fichiers MiSeq: HS24-BC2_S33_L001_R1_001.fastq.gz -> HS24-BC2
    # Pour les fichiers iSeq: HS24-BC2_R1.fastq -> HS24-BC2
    if [[ "$filename" =~ ^([^_]+)_(S[0-9]+_L[0-9]+_)?R[12](_[0-9]+)?\.(fastq|fq)(\.gz)?$ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "UNKNOWN"
    fi
}

# Function pour détecter le type de séquenceur
detect_sequencer_type() {
    local filename="$1"
    if [[ "$filename" =~ _S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz$ ]]; then
        echo "MiSeq"
    elif [[ "$filename" =~ _R[12]\.fastq$ ]]; then
        echo "iSeq"
    else
        echo "Unknown"
    fi
}

# Function pour détecter si c'est un contrôle
is_control_sample() {
    local sample_name="$1"
    # Patterns pour détecter les contrôles (case insensitive)
    if [[ "$sample_name" =~ ^(Eau|Water|Blank|Control|Neg|Negative|NC|NTC)$ ]]; then
        return 0  # C'est un contrôle
    else
        return 1  # Ce n'est pas un contrôle
    fi
}

log "Génération automatique des fichiers manifest et metadata"
log "Analyse du répertoire: $RAW_DATA_DIR"

cd "$RAW_DATA_DIR"

# Arrays pour stocker les informations
declare -A samples_main
declare -A samples_control
declare -A sample_sequencer
declare -A sample_paths_r1
declare -A sample_paths_r2

# Scan des fichiers R1
for r1_file in *_R1*.fastq* *_R1*.fq*; do
    [[ -f "$r1_file" ]] || continue
    
    # Trouver le fichier R2 correspondant
    r2_file="${r1_file/_R1/_R2}"
    [[ -f "$r2_file" ]] || continue
    
    sample_name=$(extract_sample_name "$r1_file")
    sequencer=$(detect_sequencer_type "$r1_file")
    
    abs_path_r1="$RAW_DATA_DIR/$r1_file"
    abs_path_r2="$RAW_DATA_DIR/$r2_file"
    
    log "Trouvé: $sample_name ($sequencer) - $r1_file / $r2_file"
    
    # Stocker les informations
    sample_sequencer["$sample_name"]="$sequencer"
    sample_paths_r1["$sample_name"]="$abs_path_r1"
    sample_paths_r2["$sample_name"]="$abs_path_r2"
    
    # Classer en échantillon principal ou contrôle
    if is_control_sample "$sample_name"; then
        samples_control["$sample_name"]=1
        log "  -> Classé comme CONTRÔLE"
    else
        samples_main["$sample_name"]=1
        log "  -> Classé comme ÉCHANTILLON"
    fi
done

# Génération du manifest principal
log "Génération du manifest principal"
manifest_file="${DATABASE_DIR}/manifest"
cat > "$manifest_file" << EOF
sample-id	forward-absolute-filepath	reverse-absolute-filepath
EOF

for sample in "${!samples_main[@]}"; do
    echo -e "${sample}\t${sample_paths_r1[$sample]}\t${sample_paths_r2[$sample]}" >> "$manifest_file"
done

log "Manifest principal créé: $manifest_file ($(( ${#samples_main[@]} )) échantillons)"

# Génération du manifest contrôle (si des contrôles existent)
if [ ${#samples_control[@]} -gt 0 ]; then
    log "Génération du manifest contrôle"
    manifest_control_file="${DATABASE_DIR}/manifest_control"
    cat > "$manifest_control_file" << EOF
sample-id	forward-absolute-filepath	reverse-absolute-filepath
EOF
    
    for sample in "${!samples_control[@]}"; do
        echo -e "${sample}\t${sample_paths_r1[$sample]}\t${sample_paths_r2[$sample]}" >> "$manifest_control_file"
    done
    
    log "Manifest contrôle créé: $manifest_control_file ($(( ${#samples_control[@]} )) contrôles)"
else
    log "Aucun contrôle détecté - pas de manifest_control généré"
fi

# Génération des métadonnées
log "Génération du fichier sample-metadata.tsv"
metadata_file="${DATABASE_DIR}/sample-metadata.tsv"
cat > "$metadata_file" << EOF
sample-id	sequencer	sample-type	description
EOF

# Ajouter les échantillons principaux
for sample in "${!samples_main[@]}"; do
    sequencer="${sample_sequencer[$sample]}"
    echo -e "${sample}\t${sequencer}\tsample\tEnvironmental sample" >> "$metadata_file"
done

# Ajouter les contrôles
for sample in "${!samples_control[@]}"; do
    sequencer="${sample_sequencer[$sample]}"
    echo -e "${sample}\t${sequencer}\tcontrol\tNegative control" >> "$metadata_file"
done

log "Metadata créé: $metadata_file"

# Résumé
log "=== RÉSUMÉ ==="
log "Échantillons principaux: ${#samples_main[@]}"
for sample in "${!samples_main[@]}"; do
    log "  - $sample (${sample_sequencer[$sample]})"
done

if [ ${#samples_control[@]} -gt 0 ]; then
    log "Contrôles: ${#samples_control[@]}"
    for sample in "${!samples_control[@]}"; do
        log "  - $sample (${sample_sequencer[$sample]})"
    done
fi

log "Fichiers générés:"
log "  - $manifest_file"
[ ${#samples_control[@]} -gt 0 ] && log "  - $manifest_control_file"
log "  - $metadata_file"

log "Génération terminée avec succès!"