#!/bin/bash 
#SBATCH --account=OD-229285
#SBATCH --job-name build_kraken
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 24:00:00
#SBATCH --partition io  

set -e
eval "$(conda shell.bash hook)"
conda activate kraken2

workingdir=$(pwd)
mkdir -p "$workingdir/krakendb"

# Step 1: Download taxonomy and RefSeq viral sequences
kraken2-build --download-taxonomy --db "$workingdir/krakendb"
echo " Download taxonomy successful"

kraken2-build --download-library viral --db "$workingdir/krakendb"
echo " Download viral RefSeq successful"

# Step 2: Extract unique species names from existing RefSeq viral metadata
refseq_summary="$workingdir/krakendb/library/viral/assembly_summary.txt"
refseq_taxid_file="$workingdir/krakendb/taxonomy/refseq_viral_taxids.txt"

if [[ -f "$refseq_summary" ]]; then
  echo " Extracting TaxIDs from RefSeq metadata..."
  cut -f6 "$refseq_summary" | tail -n +2 | sort | uniq > "$refseq_taxid_file"
  echo " Saved TaxID exclusion list to $refseq_taxid_file"
else
  echo " RefSeq metadata file not found: $refseq_summary"
  exit 1
fi

# Step 3: Run Python downloader
cp download_viral_spp_nuccore.py "$workingdir/krakendb"
cd "$workingdir/krakendb"
python download_viral_spp_nuccore.py
echo " Python sequence downloader run successfully"

# Step 4: Final Kraken2 DB build
kraken2-build --add-to-library viral_combined_sequences_kraken.fasta --db "$workingdir/krakendb"
echo " Added sequences to Kraken2 library"

kraken2-build --build --db "$workingdir/krakendb" --threads 1 --kmer-len 30 --minimizer-len 29 --minimizer-spaces 7
