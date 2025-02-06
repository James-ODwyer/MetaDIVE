#!/bin/bash 
#SBATCH --account=XXXXXXX              # Required for CSIRO HPC. You need to specify your account. To see yours, write into a terminal (putty or other) get_project_codes
#SBATCH --job-name SampleX      # named whatever you would like but I usually name it related to what I'm analysing 
#SBATCH --nodes 1                        # nodes to use on the hpc, 1 node= max 64CPUs so leave as 1. 
#SBATCH --ntasks-per-node 1              # ntasks per node (not needed to play around with for the pipeline. the pipeline will allocate all resources as best needed)
#SBATCH --cpus-per-task 48               # total number of CPUs to allocate. depending on size of data and urgency, 12-48  
#SBATCH --mem 100G                       # Total memory. Can require a lot particularly if you want to run trinity! between 80 and 180 depending on complexity of data
#SBATCH --time 48:00:00                 # Time requirements hh/mm/ss would recommend around 100 hours for large datasets. if it doesn't complete you can always launch the script again

# Need to activate conda through source when running a slurm script

eval "$(conda shell.bash hook)"
#conda activate worked from snakemake7 preloaded in interactive. Will need to update this with your conda path! Same thing as below. Update the conda prefix! 
conda activate snakemake7

# Quickly unlock working directory in case its locked (can happen if you launch and then cancel a run mid way through) 
snakemake -s snakefile.snakefile -j ${SLURM_CPUS_PER_TASK} --unlock --use-conda
#run the main pipeline
snakemake -s snakefile.snakefile -j ${SLURM_CPUS_PER_TASK} --keep-going --use-conda --rerun-triggers mtime --rerun-incomplete -T 3 --resources genbank=3 trinity=2


# this is how to get the full summary results. I haven't incorporated this in fully yet so after the run has completed you will need to run these two lines in the terminal with your working directory 

#cp scripts/extract_results.sh .
#bash extract_results.sh