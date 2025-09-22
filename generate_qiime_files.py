#!/usr/bin/env python3
"""
Script automatisé pour générer les fichiers manifest et metadata QIIME2
Usage: python3 generate_qiime_files.py
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Configuration
ROOTDIR = "/nvme/bio/data_fungi/valormicro_nc"
RAW_DATA_DIR = f"{ROOTDIR}/01_raw_data"
DATABASE_DIR = f"{ROOTDIR}/98_databasefiles"

# Patterns pour détecter les contrôles (insensible à la casse)
CONTROL_PATTERNS = [
    r'^eau$', r'^water$', r'^blank$', r'^control$', 
    r'^neg$', r'^negative$', r'^nc$', r'^ntc$',
    r'^ctrl', r'^blank[0-9]*$', r'^neg[0-9]*$'
]

def extract_sample_name(filename):
    """Extrait le nom d'échantillon depuis le nom de fichier"""
    # Pour MiSeq: HS24-BC2_S33_L001_R1_001.fastq.gz -> HS24-BC2
    # Pour iSeq: HS24-BC2_R1.fastq -> HS24-BC2
    # Gérer aussi les cas spéciaux comme POE2 vs Poe2MN
    
    patterns = [
        r'^([^_]+)_S\d+_L\d+_R[12]_\d+\.(fastq|fq)(\.gz)?$',  # MiSeq
        r'^([^_]+)_R[12]\.(fastq|fq)(\.gz)?$',                # iSeq simple
        r'^([^_]+)_.*_R[12]\.(fastq|fq)(\.gz)?$'              # Autres formats
    ]
    
    for pattern in patterns:
        match = re.match(pattern, filename, re.IGNORECASE)
        if match:
            sample_name = match.group(1)
            # Normaliser les noms similaires
            sample_name = sample_name.replace('Poe2MN', 'POE2')
            return sample_name
    
    return "UNKNOWN"

def detect_sequencer_type(filename):
    """Détecte le type de séquenceur basé sur le format du nom"""
    if re.search(r'_S\d+_L\d+_R[12]_\d+\.fastq\.gz$', filename, re.IGNORECASE):
        return "MiSeq"
    elif re.search(r'_R[12]\.fastq$', filename, re.IGNORECASE):
        return "iSeq"
    else:
        return "Unknown"

def is_control_sample(sample_name):
    """Détermine si l'échantillon est un contrôle"""
    for pattern in CONTROL_PATTERNS:
        if re.match(pattern, sample_name, re.IGNORECASE):
            return True
    return False

def find_fastq_pairs(raw_data_dir):
    """Trouve tous les paires de fichiers FASTQ R1/R2"""
    pairs = {}
    
    # Chercher tous les fichiers R1
    for filepath in Path(raw_data_dir).glob("*_R1*"):
        if not filepath.is_file():
            continue
            
        filename = filepath.name
        # Construire le nom du fichier R2 correspondant
        r2_filename = filename.replace('_R1', '_R2')
        r2_filepath = filepath.parent / r2_filename
        
        if r2_filepath.exists():
            sample_name = extract_sample_name(filename)
            sequencer = detect_sequencer_type(filename)
            
            pairs[sample_name] = {
                'r1_path': str(filepath),
                'r2_path': str(r2_filepath),
                'sequencer': sequencer,
                'is_control': is_control_sample(sample_name)
            }
            
            print(f"Trouvé: {sample_name} ({sequencer}) - {filename} / {r2_filename}")
            if pairs[sample_name]['is_control']:
                print(f"  -> Classé comme CONTRÔLE")
            else:
                print(f"  -> Classé comme ÉCHANTILLON")
        else:
            print(f"ATTENTION: Fichier R2 manquant pour {filename}")
    
    return pairs

def generate_manifest(samples, output_path, include_controls=False):
    """Génère un fichier manifest QIIME2"""
    with open(output_path, 'w') as f:
        f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        
        count = 0
        for sample_name, info in samples.items():
            if include_controls == info['is_control']:
                f.write(f"{sample_name}\t{info['r1_path']}\t{info['r2_path']}\n")
                count += 1
        
        print(f"Manifest créé: {output_path} ({count} échantillons)")

def generate_metadata(samples, output_path):
    """Génère un fichier metadata QIIME2"""
    data = []
    
    for sample_name, info in samples.items():
        sample_type = "control" if info['is_control'] else "sample"
        description = "Negative control" if info['is_control'] else "Environmental sample"
        
        data.append({
            'sample-id': sample_name,
            'sequencer': info['sequencer'],
            'sample-type': sample_type,
            'description': description
        })
    
    # Créer DataFrame et sauvegarder
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Metadata créé: {output_path}")

def main():
    """Fonction principale"""
    print("=== Générateur automatique de fichiers QIIME2 ===")
    print(f"Analyse du répertoire: {RAW_DATA_DIR}")
    
    # Créer le répertoire de base de données
    Path(DATABASE_DIR).mkdir(parents=True, exist_ok=True)
    
    # Trouver tous les pairs FASTQ
    samples = find_fastq_pairs(RAW_DATA_DIR)
    
    if not samples:
        print("ERREUR: Aucun fichier FASTQ trouvé!")
        sys.exit(1)
    
    # Séparer les échantillons et contrôles
    main_samples = {k: v for k, v in samples.items() if not v['is_control']}
    control_samples = {k: v for k, v in samples.items() if v['is_control']}
    
    # Générer le manifest principal
    manifest_path = f"{DATABASE_DIR}/manifest"
    generate_manifest(main_samples, manifest_path, include_controls=False)
    
    # Générer le manifest contrôle si nécessaire
    if control_samples:
        manifest_control_path = f"{DATABASE_DIR}/manifest_control"
        generate_manifest(control_samples, manifest_control_path, include_controls=True)
    else:
        print("Aucun contrôle détecté - pas de manifest_control généré")
    
    # Générer les métadonnées
    metadata_path = f"{DATABASE_DIR}/sample-metadata.tsv"
    generate_metadata(samples, metadata_path)
    
    # Résumé
    print("\n=== RÉSUMÉ ===")
    print(f"Échantillons principaux: {len(main_samples)}")
    for sample in main_samples.keys():
        print(f"  - {sample} ({samples[sample]['sequencer']})")
    
    if control_samples:
        print(f"Contrôles: {len(control_samples)}")
        for sample in control_samples.keys():
            print(f"  - {sample} ({samples[sample]['sequencer']})")
    
    print("\nFichiers générés:")
    print(f"  - {manifest_path}")
    if control_samples:
        print(f"  - {DATABASE_DIR}/manifest_control")
    print(f"  - {metadata_path}")
    
    print("\nGénération terminée avec succès!")

if __name__ == "__main__":
    main()