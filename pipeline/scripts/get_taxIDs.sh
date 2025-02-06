#!/bin/bash 


while getopts 'a:i:o:p:' c
do
  case $c in
    a) diamond_output=$OPTARG ;;
    i) taxonkit_db=$OPTARG ;;
    o) updated_output=$OPTARG ;;
    p) threads=$OPTARG ;;
  esac
done



# Script to add missing taxids into correct spots using TaxonKit
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo " starting first stage of identification $current_date_time "



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


# Step 3: Parallelized TaxonKit search
export taxonkit_db

process_species() {
    local species="$1"
    new_taxid=$(echo "$species" | taxonkit name2taxid --data-dir "$taxonkit_db" | cut -f2)
    if [ -n "$new_taxid" ]; then
        echo -e "$species\t$new_taxid"
    fi
}

export -f process_species

parallel --jobs $threads process_species {} "$taxonkit_db" :::: "$species_file.uniq" >> "$taxid_map_file"

# Step 4: Create dictionary array from taxid_map_file
declare -A taxid_map
while IFS=$'\t' read -r species taxid; do
    taxid_map["$species"]="$taxid"
done < "$taxid_map_file"

# Step 5: Apply TaxonKit results to relevant rows
while IFS=$'\t' read -r species old_taxid row_num rest_of_line; do
    IFS=$'\t' read -r -a fields <<< "$rest_of_line"
    if [ -n "${taxid_map["$species"]}" ]; then
        taxid="${taxid_map["$species"]}"
        fields+=("$taxid")
        (IFS=$'\t'; echo "${fields[*]}") >> "$temp_output"
        ((new_taxid_count++))
    fi
done < "$species_file"

# Handle rows that already had valid TaxID initially
while IFS=$'\t' read -r -a line; do
    taxid="${line[6]}"
    if [[ "$taxid" =~ ^([0-9]+(:[0-9]+)*)(;[0-9]+(:[0-9]+)*)*$ ]]; then
        line+=("$taxid")
        (IFS=$'\t'; echo "${line[*]}") >> "$temp_output"
    fi
done < "$processed_output"

# Move the temporary file to the final output
mv "$temp_output" "$updated_output"

# Clean up temporary files
rm "$species_file" "$species_file.uniq" "$taxid_map_file" "$processed_output"

# Print the count of new TaxIDs added
echo "Number of new TaxIDs added: $new_taxid_count"
echo "Updated Diamond output written to $updated_output"

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo "Finished first stage of identification: $current_date_time"






















