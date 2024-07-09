#!/bin/bash 
#SBATCH --account=OD-229285
#SBATCH --job-name TaxonKit_IDs
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 1:00:00



while getopts 'a:i:o:' c
do
  case $c in
    a) diamond_output=$OPTARG ;;
    i) taxonkit_db=$OPTARG ;;
    o) updated_output=$OPTARG ;;
  esac
done



# Script to add missing taxids into correct spots using TaxonKit


# Note, still need to confirm argparse is downloaded in RdataplottingS3.
# Note,  will need to touch an output file at the end of the script and generate a while wait loop inside the snakemake rule.
# Define file paths (Note, still need to incorporate the snakemake variables to be imported into the script)
# Need 1. The diamond output file
# 2. The filepath to the TaxonKit database (just the names and nodes.dmp files extracted directly from NCBI using 
#wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
#tar -zxvf taxdump.tar.gz
# 3. The output file (note Will need to replace this with the snakemake output file names)

# Components of the Diamond output that are important 
# column 2=sseqid code 
# column 7= taxid code
# column8=stitle code

# Temporary files for holding dictionary values plus storing updated rows before moving to output file name.
temp_output=$(mktemp)
species_file=$(mktemp)
taxid_map_file=$(mktemp)
processed_output=$(mktemp)

# Step 0: Preprocess the Diamond output to replace double tabs with tab-NA-tab
# Step 0 is required because Diamond outputs empty for the taxids it doesn't find e.g., it results in \t\t instead of \t${taxid}\t
sed 's/\t\t/\tNA\t/g' "$diamond_output" > "$processed_output"

# Initialize counter for new TaxIDs. 
new_taxid_count=0

# Extract species name from stitle of Diamond file
# Name always located within the square brackets, easy to locate and extract using regex
# Create function to do so repeatedly. 
extract_species_name() {
    local title="$1"
    echo "$title" | grep -oP "\[\K[^\]]+"
}

# Step 1: Identify rows missing taxonomy data and store their row numbers
i=0
while IFS=$'\t' read -r -a line; do
    accession="${line[1]}"
    taxid="${line[6]}"
    title="${line[7]}"
    # taxID's can take the form of one or multiple strings of numbers separated by ';' hence the requirement for a more complex regex pattern.
    # Read all tables in as tab and export them out as tab delimited (important because whitespace exists within $stitle
    if [[ ! "$taxid" =~ ^([0-9]+(:[0-9]+)*)(;[0-9]+(:[0-9]+)*)*$ ]]; then
        species=$(extract_species_name "$title")
        if [ -n "$species" ]; then
            (IFS=$'\t'; echo -e "$species\t$taxid\t$i\t${line[*]}") >> "$species_file"
        fi
    fi
    ((i++))
done < "$processed_output"

# Step 2: Extract and collect unique species names
cut -f1 "$species_file" | sort | uniq > "$species_file.uniq"

# Step 3: Search NCBI database (through TaxonKit) once for each unique species name
# TaxonKit appears to be quite slow (~3 seconds per index search)
# Taxonomizer in R is <1 second per search marking a significant speed boost likely due to the optimised sql indexing of the reference database. So will just determine taxID's in this script and leave wider taxonomy to importing into R.
while read -r species; do
    new_taxid=$(echo "$species" | taxonkit name2taxid --data-dir "$taxonkit_db" | cut -f2)
    if [ -n "$new_taxid" ]; then
        (IFS=$'\t'; echo -e "$species\t$new_taxid") >> "$taxid_map_file"
    fi
done < "$species_file.uniq"

# Create dictionary array from taxid_map_file
declare -A taxid_map
while IFS=$'\t' read -r species taxid; do
    taxid_map["$species"]="$taxid"
done < "$taxid_map_file"

# Step 4: Apply TaxonKit results to relevant rows
while IFS=$'\t' read -r species old_taxid row_num rest_of_line; do
    IFS=$'\t' read -r -a fields <<< "$rest_of_line"
    if [ -n "${taxid_map["$species"]}" ]; then
        taxid="${taxid_map["$species"]}"
        fields+=("$taxid")
        # Reconstruct the line with tabs
        (IFS=$'\t'; echo "${fields[*]}") >> "$temp_output"
        ((new_taxid_count++))
#    else
#        fields+=("NoTaxIDFound") 
#        # Reconstruct the line with tabs
#        (IFS=$'\t'; echo "${fields[*]}") >> "$temp_output"
    fi
done < "$species_file"
# Need to ignore the rows already identified above. 
# Handle rows that had valid TaxID initially
while IFS=$'\t' read -r -a line; do
    taxid="${line[6]}"
    if [[ "$taxid" =~ ^([0-9]+(:[0-9]+)*)(;[0-9]+(:[0-9]+)*)*$ ]]; then
        line+=("$taxid")
        (IFS=$'\t'; echo "${line[*]}") >> "$temp_output"
#    else
#        line+=("NoTaxIDFound")
    fi
    # Reconstruct the line with tabs

done < "$processed_output"

# Move the temporary file to the final output
mv "$temp_output" "$updated_output"

# Clean up temporary files
rm "$species_file" "$species_file.uniq" "$taxid_map_file" "$processed_output"

# Print the count of new TaxIDs added
echo "Number of new TaxIDs added: $new_taxid_count"
echo "Updated Diamond output written to $updated_output"



























