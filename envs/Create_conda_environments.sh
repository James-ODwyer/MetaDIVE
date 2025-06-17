#!/bin/bash
#SBATCH --account=OD-229285
#SBATCH --job-name building_conda_environments     # named whatever you would like
#SBATCH --nodes 1                        # nodes to use on the hpc
#SBATCH --ntasks-per-node 1              # ntasks per node (not needed to play around with for the pipeline. the pipeline will allocate all resources as best needed)
#SBATCH --cpus-per-task 1               # total number of CPUs to allocate
#SBATCH --mem 4G                       # Total memory.
#SBATCH --time 2:00:00                 # Time requirements hh/mm/ss
#SBATCH --partition io                # If a particular download node is required


# Internet access is needed

cd ./envs

conda_path=$(conda info --base)
conda_base=$(conda info --base)
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

mkdir -p "$conda_path"/envs/kraken2

env_path="$conda_path"/envs/kraken2

# Create the Conda environment from the YAML file at the specified location
conda env create --file kraken2.yaml --prefix $env_path

echo "Conda environment created at: $env_path"
echo "Conda environment kraken2 fully installed"


# Collect all .condarc paths from config show-sources
condarc_files=$(conda config --show-sources | awk '/^==> / {print $2}')

# Loop over all discovered .condarc files
for file in $condarc_files; do
    echo " Updating .condarc: $file"

    # Make sure the file exists (it should, but for safety)
    [ -f "$file" ] || touch "$file"

    # Patch envs_dirs
    if ! grep -q "^envs_dirs:" "$file"; then
        echo -e "\nenvs_dirs:\n  - ${conda_base}\n  - ${conda_base}/envs" >> "$file"
        echo "  Added envs_dirs to $file"
    else
        grep -q "  - ${conda_base}" "$file" || sed -i "/^envs_dirs:/a \  - ${conda_base}" "$file"
        grep -q "  - ${conda_base}/envs" "$file" || sed -i "/^envs_dirs:/a \  - ${conda_base}/envs" "$file"
    fi

    # Patch pkgs_dirs
    if ! grep -q "^pkgs_dirs:" "$file"; then
        echo -e "\npkgs_dirs:\n  - ${conda_base}/pkgs" >> "$file"
        echo "  Added pkgs_dirs to $file"
    else
        grep -q "  - ${conda_base}/pkgs" "$file" || sed -i "/^pkgs_dirs:/a \  - ${conda_base}/pkgs" "$file"
    fi
done

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

eval "$(conda shell.bash hook)"

conda activate Rdataplotting

Rscript install_d3Tree.R

conda deactivate
