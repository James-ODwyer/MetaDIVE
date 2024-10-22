#!/bin/bash 
#SBATCH --account=OD-229285              # Required for CSIRO HPC. You need to specify your account. To see yours, write into a terminal (putty or other) get_project_codes
#SBATCH --job-name build_kraken      # named whatever you would like but I usually name it related to what I'm analysing 
#SBATCH --nodes 1                        # nodes to use on the hpc, 1 node= max 64CPUs so leave as 1. 
#SBATCH --ntasks-per-node 1              # ntasks per node (not needed to play around with for the pipeline. the pipeline will allocate all resources as best needed)
#SBATCH --cpus-per-task 2               # total number of CPUs to allocate. depending on size of data and urgency, 12-48  
#SBATCH --mem 12G                       # Total memory. Can require a lot particularly if you want to run trinity! between 80 and 180 depending on complexity of data
#SBATCH --time 4:00:00                 # Time requirements hh/mm/ss would recommend around 100 hours for large datasets. if it doesn't complete you can always launch the script again
#SBATCH --partition io  



eval "$(conda shell.bash hook)"

conda activate kraken2


workingdir=$(pwd)

mkdir "$workingdir"/krakendb

kraken2-build --download-taxonomy --db "$workingdir"/krakendb

echo " Download taxonomy run succesfully"

kraken2-build --download-library viral --db "$workingdir"/krakendb

echo " Download viral refseq run succesfully"

cp download_viral_spp_nuccore.py "$workingdir"/krakendb
cd "$workingdir"/krakendb

python "$workingdir"/download_viral_spp_nuccore.py

echo " python run succesfully"

kraken2-build --add-to-library viral_sequences_selected.fasta --db "$workingdir"/krakendb

echo " add to library run succesfully"
kraken2-build --build --db "$workingdir"/krakendb --threads 2 --kmer-len 29 --minimizer-len 29 --minimizer-spaces 7
