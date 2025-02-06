#!/bin/bash 
#SBATCH --account=OD-229285
#SBATCH --job-name Compile_results_pipeline
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 2:00:00


# Summary results extraction script. 
# This will not work effectively if you change the names of program directories (e.g., sub directory renames) 
programdir=$(pwd)

mkdir "$programdir"/Summarised_results
mkdir "$programdir"/Summarised_results/CO1_results
mkdir "$programdir"/Summarised_results/CO1_results/assembled_contigs
mkdir "$programdir"/Summarised_results/CO1_results/blastn_contigs
mkdir "$programdir"/Summarised_results/CO1_results/read_homologies
mkdir "$programdir"/Summarised_results/LSU_results
mkdir "$programdir"/Summarised_results/LSU_results/assembled_contigs
mkdir "$programdir"/Summarised_results/LSU_results/blastn_contigs
mkdir "$programdir"/Summarised_results/LSU_results/read_homologies
mkdir "$programdir"/Summarised_results/SSU_results
mkdir "$programdir"/Summarised_results/SSU_results/assembled_contigs
mkdir "$programdir"/Summarised_results/SSU_results/blastn_contigs
mkdir "$programdir"/Summarised_results/SSU_results/read_homologies
mkdir "$programdir"/Summarised_results/Reference_guided_viral_genomes
mkdir "$programdir"/Summarised_results/Reference_guided_viral_genomes/assembled_genomes_and_reads
mkdir "$programdir"/Summarised_results/Reference_guided_viral_genomes/mafft_alignments
mkdir "$programdir"/Summarised_results/Reference_guided_viral_genomes/finished_trees
mkdir "$programdir"/Summarised_results/Summary_results_per_sample_contigs
mkdir "$programdir"/Summarised_results/Summary_results_combined_figures_and_tables_contigs
mkdir "$programdir"/Summarised_results/raw_assemblies
mkdir "$programdir"/Summarised_results/Metabat_binned_contigs
mkdir "$programdir"/Summarised_results/combined_raw_reads_and_contigs_viruses


# Need to update for future pipelines as the blastn '/' error is fixed now so the results are saved to 04n sub directory

# If statement to detect if the microbiome was run? 

microbiomepresence="$programdir/03_CO1_FILTER/HOSTREMOVAL/03_LCA_RESULTS/*host_species*"
if ls $microbiomepresence 1> /dev/null 2>&1; then

# Do CO1
cp "$programdir"/03_CO1_FILTER/HOSTREMOVAL/03_LCA_RESULTS/*.html "$programdir"/Summarised_results/CO1_results/read_homologies
# Copy and rename the contigs files
source_path_pattern="$programdir/03_CO1_FILTER/HOSTREMOVAL/01a_STRONG_ALIGN_CONTIGS/*/final.contigs.fa"
for file in $source_path_pattern; do
  # Extract the sample name from the path
  sample_name=$(basename $(dirname "$file"))
  
  # Define the new file name
  new_file_name="${sample_name}_final.contigs.fa"
  
  # Copy and rename the file
  cp "$file" "$programdir"/Summarised_results/CO1_results/blastn_contigs/"$new_file_name"
  
done

cp "$programdir"/03_CO1_FILTER/HOSTREMOVAL/04_CONTIG_BLASTN_RESULTS/*Contigsall* "$programdir"/Summarised_results/CO1_results/assembled_contigs


# Do LSU
cp "$programdir"/04_LSU_FILTER/HOSTREMOVAL/03_LCA_RESULTS/*.html "$programdir"/Summarised_results/LSU_results/read_homologies

# Copy and rename the contigs files
source_path_pattern="$programdir/04_LSU_FILTER/HOSTREMOVAL/01a_STRONG_ALIGN_CONTIGS/*/final.contigs.fa"
for file in $source_path_pattern; do
  # Extract the sample name from the path
  sample_name=$(basename $(dirname "$file"))
  
  # Define the new file name
  new_file_name="${sample_name}_final.contigs.fa"
  
  # Copy and rename the file
  cp "$file" "$programdir"/Summarised_results/LSU_results/blastn_contigs/"$new_file_name"
done

cp "$programdir"/05_SSU_FILTER/HOSTREMOVAL/04_CONTIG_BLASTN_RESULTS/*Contigsall* "$programdir"/Summarised_results/SSU_results/assembled_contigs

# Do SSU
cp "$programdir"/05_SSU_FILTER/HOSTREMOVAL/03_LCA_RESULTS/*.html "$programdir"/Summarised_results/SSU_results/read_homologies

# Copy and rename the contigs files
source_path_pattern="$programdir/05_SSU_FILTER/HOSTREMOVAL/01a_STRONG_ALIGN_CONTIGS/*/final.contigs.fa"
for file in $source_path_pattern; do
  # Extract the sample name from the path
  sample_name=$(basename $(dirname "$file"))
  
  # Define the new file name
  new_file_name="${sample_name}_final.contigs.fa"
  
  # Copy and rename the file
  cp "$file" "$programdir"/Summarised_results/SSU_results/blastn_contigs/"$new_file_name"
  
done

cp "$programdir"/05_SSU_FILTER/HOSTREMOVAL/04_CONTIG_BLASTN_RESULTS/*Contigsall* "$programdir"/Summarised_results/SSU_results/assembled_contigs

# Microbiome done. 

fi

# Copy over trees if tree whole genome analysis undertaken


fingenomespresence="$programdir/18B_COMPLETE_VIRAL_GENOMES_GUIDED_ASSEMBLY/*"
if ls $fingenomespresence 1> /dev/null 2>&1; then


cp -r "$programdir"/18E_IQTREE_GENERATED_VIRAL_GENOMES/*/ "$programdir"/Summarised_results/Reference_guided_viral_genomes/finished_trees

cp -r "$programdir"/18D_TRIMAL_ALIGNMENTS_GENERATED_VIRAL_GENOMES/*/ "$programdir"/Summarised_results/Reference_guided_viral_genomes/mafft_alignments

source_dir1="${programdir}/18B_COMPLETE_VIRAL_GENOMES_GUIDED_ASSEMBLY/*/moderate_coverage_genomes"
source_dir2="${programdir}/18B_COMPLETE_VIRAL_GENOMES_GUIDED_ASSEMBLY/*/finished_genomes"
target_base_dir="${programdir}/Summarised_results/Reference_guided_viral_genomes/assembled_genomes_and_reads"

# Iterate over each sample directory
for sample_dir in ${programdir}/18B_COMPLETE_VIRAL_GENOMES_GUIDED_ASSEMBLY/*; do
  # Extract the sample name
  sample_name=$(basename "$sample_dir")
  
  # Create the target directory for the sample
  target_dir="${target_base_dir}/${sample_name}"
  mkdir -p "$target_dir"
  
  # Copy the moderate_coverage_genomes directory
  if [ -d "${sample_dir}/moderate_coverage_genomes" ]; then
    cp -r "${sample_dir}/moderate_coverage_genomes" "$target_dir"
  fi
  
  # Copy the finished_genomes directory
  if [ -d "${sample_dir}/finished_genomes" ]; then
    cp -r "${sample_dir}/finished_genomes" "$target_dir"
  fi
done


fi


# summary per sample results copy
cp "$programdir"/99_SUMMARY_RESULTS/*.txt "$programdir"/Summarised_results/Summary_results_per_sample_contigs

# Don't want the raw diamonds table (they are a few GB on big runs and superfluous stats anyway 
find "$programdir/Summarised_results" -type f -name '*_raw_diamond_hits.txt' -exec rm -f {} \;

# summary all samples results copy 
cp "$programdir"/99_SUMMARY_RESULTS_COMBINED_GRAPHS/* "$programdir"/Summarised_results/Summary_results_combined_figures_and_tables_contigs

# raw assemblies for initial megahit/trinity assembly
cp "$programdir"/06_CONTIG_ASSEMBLY/* "$programdir"/Summarised_results/raw_assemblies



Metabatresults="$programdir/26_METABAT_BINS_SORTED/*"
if ls $Metabatresults1> /dev/null 2>&1; then

cp -r "$programdir"/26_METABAT_BINS_SORTED/*/ "$programdir"/Summarised_results/Metabat_binned_contigs
fi

raw_contigs_sort="$programdir/99_SUMMARY_RESULTS_RAWS_ADDED_TABLES/*"
if ls $raw_contigs_sort> /dev/null 2>&1; then

cp -r "$programdir"/99_SUMMARY_RESULTS_RAWS_ADDED_TABLES/*/ "$programdir"/Summarised_results/combined_raw_reads_and_contigs_viruses

fi


# Path to the Snakemake config file
CONFIG_FILE="config.yaml"

# Extract Delete_inter_files variable from the config
DELETE_INTER_FILES=$(grep -oP '(?<=Delete_inter_files: ")[^"]*' "$CONFIG_FILE" || echo "no")

# Conditional statement to check value
if [[ "$DELETE_INTER_FILES" == "yes" ]]; then
    echo "Deleting intermediary files..."
    #rm -rf ./path/to/intermediary/files/*
elif [[ "$DELETE_INTER_FILES" == "no" ]]; then
    echo "No intermediary files will be deleted."
else
    echo "No intermediary files will be deleted (invalid value in config)."
fi


