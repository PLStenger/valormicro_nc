# valormicro_nc
Characterization of marine microbial resources for analysis and enhancement of New Caledonia's natural heritage - Project from Drs **Anton Véronique** (CNRS/IRD Nouméa - New Caledonia)

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/valormicro_nc.git

### Update the pipeline in local by :

    git pull
    
# VALORMICRO PIPELINE INSTALLATION - DIRECT COMMANDS

# 1. Create and navigate to working directory
export ROOTDIR="/XXXXXX/valormicro_nc"
mkdir -p "$ROOTDIR"/{00_scripts,01_raw_data,02_qualitycheck,03_cleaned_data,04_qualitycheck,05_QIIME2/{core,visual},98_databasefiles,99_softwares/adapters,exports,tmp}
cd "$ROOTDIR"

# 2. Install Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p "$HOME/miniconda3"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda init bash
source ~/.bashrc

# 3. Install Mamba to accelerate installations
conda install -y -c conda-forge mamba

# 4. Create conda environments
mamba create -y -n fastqc -c bioconda fastqc
mamba create -y -n multiqc -c bioconda multiqc  
mamba create -y -n trimmomatic -c bioconda trimmomatic
mamba create -y -n python-bio python=3.9 pandas numpy

# 5. Install QIIME2 2021.4
wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-linux-conda.yml
mamba env create -n qiime2-2021.4 --file qiime2-2021.4-py38-linux-conda.yml
rm qiime2-2021.4-py38-linux-conda.yml

# 6. Install GNU Parallel
mamba install -y -c conda-forge parallel

# 7. Create Trimmomatic adapter file
cat > "$ROOTDIR/99_softwares/adapters/sequences.fasta" << 'EOF'
>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_f
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_PE_r
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>Nextera_LongInsert_R1
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Nextera_LongInsert_R2
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
EOF

# 8. Download SILVA database
cd "$ROOTDIR/98_databasefiles"
wget -O silva-138-99-515-926-nb-classifier.qza \
    "https://data.qiime2.org/2021.4/common/silva-138-99-515-806-nb-classifier.qza"

# 9. Add configuration to .bashrc
cat >> ~/.bashrc << 'EOF'

# === VALORMICRO PIPELINE ===
export VALORMICRO_ROOT="/nvme/bio/data_fungi/valormicro_nc"
alias qiime2='conda activate qiime2-2021.4'
alias fastqc='conda activate fastqc'
EOF

source ~/.bashrc




# Test environments
conda env list
conda activate qiime2-2021.4 && qiime --version
conda activate fastqc && fastqc --version  
conda activate trimmomatic && trimmomatic -version
parallel --version

# Test files
ls -la "$ROOTDIR/99_softwares/adapters/"
ls -la "$ROOTDIR/98_databasefiles/"



# If conda is not found
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# If permission issues
sudo chown -R $(whoami):$(whoami) "$ROOTDIR"

# If QIIME2 fails to install
mamba create -n qiime2-2021.4 -c qiime2 -c bioconda -c conda-forge qiime2=2021.4

# If SILVA database fails
wget -O "$ROOTDIR/98_databasefiles/silva-138-99-515-926-nb-classifier.qza" \
    "https://data.qiime2.org/2021.4/common/silva-138-99-nb-classifier.qza"
