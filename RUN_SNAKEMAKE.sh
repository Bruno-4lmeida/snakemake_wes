### CHECK LIST:
# 1 Config file and paths, make sure you review all the paths and reference genome (check if numeric)
# 2 Never edit the original, make sure to copy a config_file into confifiles dir, check this version with a --dru-run
# 3 Always check the git branch before run or edit anything. Always edit the branch dev: git switch dev


### TO RUN SNAKEMAKE GENOME WORKFLOW

cd Interomics/workflows/interomics_snakemake
conda activate interomics_snakemake

# Dry run to check config and paths
snakemake --cores all --use-conda --conda-frontend conda --keep-going --dry-run
### CHECK LIST:
# 1 Config file and paths, make sure you review all the paths and reference genome (check if numeric)
# 2 Never edit the original, make sure to copy a config_file into confifiles dir, check this version with a --dru-run
# 3 Always check the git branch before run or edit anything. Always edit the branch dev: git switch dev
# 4 Remember to add the name of your NEW config to each code below


### TO RUN SNAKEMAKE GENOME WORKFLOW
# DIR, ENV AND GIT BRANCH (use this)
cd Interomics/workflows/interomics_snakemake
conda activate interomics_snakemake
git init
git switch dev

# Dry run to check config and paths
snakemake --cores all --use-conda --conda-frontend conda --keep-going --dry-run --configfile configfiles/config.yaml # modify config

# First run
snakemake --cores all --use-conda --conda-frontend conda --keep-going --configfile configfiles/config.yaml # modify config

# Rerun
snakemake --cores all --use-conda --conda-frontend conda --keep-going --configfile configfiles/config.yaml --rerun-incomplete # modify config

# If anything was edited:
snakemake --cores all --use-conda --conda-frontend conda --keep-going --rerun-incomplete --configfile configfiles/config.yaml --touch # modify config

# First run
snakemake --cores all --use-conda --conda-frontend conda --keep-going

# Rerun
#If anything was edited:
snakemake --cores all --use-conda --conda-frontend conda --keep-going --rerun-incomplete  --touch
snakemake --cores all --use-conda --conda-frontend conda --keep-going --rerun-incomplete 

#SV
snakemake --cores all --use-conda --conda-frontend conda --snakefile Snakefile_copynumber --configfile configcnvs.yaml --keep-going --rerun-incomplete --touch