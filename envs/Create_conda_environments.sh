#!/bin/bash 
#SBATCH --account=OD-229285                         # If account info is needed
#SBATCH --job-name building_conda_environments     # named whatever you would like 
#SBATCH --nodes 1                        # nodes to use on the hpc 
#SBATCH --ntasks-per-node 1              # ntasks per node (not needed to play around with for the pipeline. the pipeline will allocate all resources as best needed)
#SBATCH --cpus-per-task 1               # total number of CPUs to allocate  
#SBATCH --mem 4G                       # Total memory. 
#SBATCH --time 6:00:00                 # Time requirements hh/mm/ss 
#SBATCH --partition io                # If a particular download node is required


# Internet access is needed

conda_path=$(conda info --base)

conda config --set channel_priority flexible

mkdir -p "$conda_path"/envs/snakemake7

env_path="$conda_path"/envs/snakemake7

# Create the Conda environment from the YAML file at the specified location
conda env create --file snakemake7.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment Snakemake7 fully installed"

echo "Updating conda.py to fix error in snakemake version 7-8"

cp conda.py "$conda_path"/envs/snakemake7/lib/python3.10/site-packages/snakemake/deployment/conda.py



mkdir -p "$conda_path"/envs/trinity

env_path="$conda_path"/envs/trinity

# Create the Conda environment from the YAML file at the specified location
conda env create --file trinity.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment trinity fully installed"




mkdir -p "$conda_path"/envs/spades

env_path="$conda_path"/envs/spades

# Create the Conda environment from the YAML file at the specified location
conda env create --file spades.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment spades fully installed"




mkdir -p "$conda_path"/envs/R4_2_dada

env_path="$conda_path"/envs/R4_2_dada

# Create the Conda environment from the YAML file at the specified location
conda env create --file R4_2_dada.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment R4_2_dada fully installed"




mkdir -p "$conda_path"/envs/Rdataplotting

env_path="$conda_path"/envs/Rdataplotting

# Create the Conda environment from the YAML file at the specified location
conda env create --file Rdataplotting.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment Rdataplotting fully installed"


mkdir -p "$conda_path"/envs/genomad

env_path="$conda_path"/envs/genomad

# Create the Conda environment from the YAML file at the specified location
conda env create --file genomad.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment genomad fully installed"

