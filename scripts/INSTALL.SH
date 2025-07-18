##### 1 INSTALL CONDA #####
# Miniconda 24.7
wget https://repo.anaconda.com/miniconda/Miniconda3-py312_24.7.1-0-Linux-x86_64.sh
bash Miniconda3-py312_24.7.1-0-Linux-x86_64.sh

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash


#### 2 CREATE SNAKEMAKE ENV ####

# Create with Mamba
conda install -n base -c conda-forge mamba
mamba env create -f /home/operator/InterOmics/workflows/interomics_snakemake/environment.yaml --channel-priority flexible

# Create with conda
conda env create -f /home/operator/InterOmics/workflows/interomics_snakemake/environment.yaml


### RUN SNAKEMAKE
conda activate interomics_snakemake

snakemake --cores all --use-conda --keep-going 