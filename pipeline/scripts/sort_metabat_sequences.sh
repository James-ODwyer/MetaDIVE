#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_metabat_sorter
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 04:00:00

while getopts 'j:p:s:z:' c
do
  case $c in
    j) genomad_results=$OPTARG ;;
    p) diamond_results=$OPTARG ;;
    s) metabat_path=$OPTARG ;;
    z) output_dir=$OPTARG ;;
  esac
done





# Iterate through Metabat2 fasta files
for metabat_fasta in "${metabat_path}"/*.fa; do
    echo "Processing ${metabat_fasta}..."
    
    # Extract Metabat2 result file name without extension
    result_file_name=$(basename "${metabat_fasta}" .fa)
    
# Initialize arrays
updated_lines=()
updated_linescontig=()
updated_linesspecies=()
updated_linesuperkingdom=()

# Iterate over each contig in the Metabat2 fasta file
while IFS= read -r contig; do
    echo "$contig"
    contig2="$(echo "$contig" | sed 's/^>//')"
    
    # Search for the contig in the Diamond results file
    hit_info=$(grep "${contig2}" "${diamond_results}" | cut -f 10)
    hit_kingdom=$(grep "${contig2}" "${diamond_results}" | cut -f 4)
    
    # Check if hit_info is not empty
    if [ -n "${hit_info}" ]; then
        updated_lines+=("${contig2} ${hit_info}")
        updated_linescontig+=("${contig2}")
        updated_linesspecies+=("${hit_info}")
        # Check if hit_kingdom is not empty
        if [ -n "${hit_kingdom}" ]; then
            echo "Adding kingdom: ${contig2} ${hit_kingdom}"
            updated_linesuperkingdom+=("${contig2} ${hit_kingdom}")
        else
            echo "Error: Kingdom information not found for ${contig2}"
        fi
    fi
done < <(grep '^>.*$' "${metabat_fasta}")

# Print results for debugging
echo "Updated Lines: ${updated_lines[@]}"
echo "Updated Line Superkingdom: ${updated_linesuperkingdom[@]}"

if [ "${#updated_lines[@]}" -gt 0 ]; then

# Still need to fix species tally. It returns the first name of the species but not full name
# Also still need to fix copying over the names

# Lastly lastly, May want to have a specific if for the unbinned or the too short because they aren't actually relevant for classifying. 
    # Tally the frequency of identified species
    species_tally=$(echo "${updated_lines[@]}" | grep -o '"[^"]*"' | sort | uniq -c | sort -nr) && \
    kingdom_tally=$(echo "${updated_linesuperkingdom[@]}" | grep -o '"[^"]*"' | sort | uniq -c | sort -nr)
    echo " kingdom tally is  $kingdom_tally"
    
        if [ "$(echo "${kingdom_tally}" | wc -l)" -eq 1 ]; then
            kingdom_name=$(echo "${kingdom_tally}" | awk '{print $2}' | tr -d '" ')
        else
            kingdom_name="mixed"
        fi
	# Take the top two most frequent species
    top_species=$(echo "${species_tally}" | cut -f 1 | grep -o '"[^"]*"' | head -n 2)

    # Concatenate the top two species names using underscores
    concatenated_species=$(echo "${top_species}" | tr '\n' '_' | tr -d '" ')

    # Create the updated contig names and sequences for the new Metabat2 results file
    updated_file_name="${output_dir}/${result_file_name}_${concatenated_species}${kingdom_name}_superkingdomkingdom.fa"

    cp ${metabat_fasta} ${updated_file_name}

for ((i=0; i<${#updated_linescontig[@]}; i++)); do
    contig="${updated_linescontig[i]}"
    species="${updated_linesspecies[i]}"
    
    # Assuming 'filename' is the name of the file to be modified
    filename="${updated_file_name}"

    # Check if the file exists before attempting modifications
    if [ -f "$filename" ]; then
        # Pattern match and replace text in the file
        sed -i "s/${contig}/&_${species}/" "$filename"
        echo "Replaced in $filename: ${contig} -> ${contig}_${species}"
    else
        echo "File not found: $filename"
    fi
done



else
    # If no hits were found, append "unassigned_contigs" to the result file name
    updated_file_name="${output_dir}/${result_file_name}_unassigned_contigs.fa"
    echo "No hits found. Writing to ${updated_file_name}..."
    cp ${metabat_fasta} ${updated_file_name}
fi

done


mkdir ${output_dir}/Eukaryota
mkdir ${output_dir}/Bacteria
mkdir ${output_dir}/Viruses
mkdir ${output_dir}/Mixed
mkdir ${output_dir}/Unassigned

mv ${output_dir}/*Eukaryota_superkingdom*fa ${output_dir}/Eukaryota
mv ${output_dir}/*Bacteria_superkingdom*fa ${output_dir}/Bacteria
mv ${output_dir}/*Viruses_superkingdom*fa ${output_dir}/Viruses
mv ${output_dir}/*mixed*fa ${output_dir}/Mixed
mv ${output_dir}/*_unassigned_*fa ${output_dir}/Unassigned




for metabat_fasta in "${output_dir}/Unassigned/*.fa"; do
    echo "Processing ${metabat_fasta}..."
    
    # Extract Metabat2 result file name without extension
    result_file_name=$(basename "${metabat_fasta}" .fa)
    
# Initialize arrays
updated_lines=()
updated_linescontig=()
updated_linestaxonomy=()

# Iterate over each contig in the Metabat2 fasta file
while IFS= read -r contig; do
    echo "$contig"
    contig2="$(echo "$contig" | sed 's/^>//')"
    
    # Search for the contig in the Diamond results file
    hit_info=$(grep "${contig2}" "${genomad_results}" | cut -f 11)
    
    # Check if hit_info is not empty
    if [ -n "${hit_info}" ]; then
        updated_lines+=("${contig2} ${hit_info}")
        updated_linescontig+=("${contig2}")
        updated_linestaxonomy+=("${hit_info}")

    fi
done < <(grep '^>.*$' "${metabat_fasta}")

filename="${metabat_fasta}"

for ((i=0; i<${#updated_linescontig[@]}; i++)); do
    contig="${updated_linescontig[i]}"
    taxonomy="${updated_linestaxonomy[i]}"
    
    # Assuming 'filename' is the name of the file to be modified
    

        # Pattern match and replace text in the file
    sed -i "s/${contig}/&_${taxonomy}/" "$filename"
    echo "Replaced in $filename: ${contig} -> ${contig}_${taxonomy}"

done

done



for metabat_fasta in "${output_dir}/Mixed/*.fa"; do
    echo "Processing ${metabat_fasta}..."
    
    # Extract Metabat2 result file name without extension
    result_file_name=$(basename "${metabat_fasta}" .fa)
    
# Initialize arrays
updated_lines=()
updated_linescontig=()
updated_linestaxonomy=()

# Iterate over each contig in the Metabat2 fasta file
while IFS= read -r contig; do
    echo "$contig"
    contig2="$(echo "$contig" | sed 's/^>//')"
    
    # Search for the contig in the Genomad results file
    hit_info=$(grep "${contig2}" "${genomad_results}" | cut -f 11)
    # Make sure that diamond results aren't being overwritten in the mixed files
    hit_infodiamond=$(grep "${contig2}" "${diamond_results}" | cut -f 10)

# Outer if statment: allows the inner if statement to be run if diamond did not identify the species
# Inner if statement: allows for the genomad results to be added into arrays to then be inserted into the fasta names. 

if [ ! -n "${hit_infodiamond}" ]; then

    # Check if hit_info is not empty
    if [ -n "${hit_info}" ]; then
        updated_lines+=("${contig2} ${hit_info}")
        updated_linescontig+=("${contig2}")
        updated_linestaxonomy+=("${hit_info}")

    fi
fi
done < <(grep '^>.*$' "${metabat_fasta}")

filename="${metabat_fasta}"

for ((i=0; i<${#updated_linescontig[@]}; i++)); do
    contig="${updated_linescontig[i]}"
    taxonomy="${updated_linestaxonomy[i]}"
    
    # Assuming 'filename' is the name of the file to be modified
    

        # Pattern match and replace text in the file
    sed -i "s/${contig}/&_${taxonomy}/" "$filename"
    echo "Replaced in $filename: ${contig} -> ${contig}_${taxonomy}"

done

done




