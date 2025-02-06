#!/bin/bash

# Parameters passed from the Snakemake rule
krakendb=$1
threads=$2
krakenreport=$3
krakenout=$4
contigs=$5
readslist=$6
log=$7
# setting additional kraken threshold to 400 contigs. Should catch almost everything 
# and not create a massive blastn backlog. Note the backlogs are about 4 hours for 3.6k contigs. therefore it should be no longer than half an hour per sample now.
reads_threshold=400

# Initialize variables to track the closest combination
closest_read_count=0
closest_conf=""
closest_hit_groups=""
closest_output=""
closest_krakenout=""
closest_krakenreport=""

# Set combinations for --minimum-hit-groups and --confidence
hit_groups=("default" "12" "18" "24")
confidences=("default" "0.5")

# Function to run Kraken2 with specific parameters
run_kraken() {
    local hit_group=$1
    local confidence=$2
    local output_suffix=$3

    # Construct Kraken2 command
    local kraken_cmd="kraken2 --db $krakendb --threads $threads \
        --report ${krakenreport%.txt}_${output_suffix}.txt \
        --output ${krakenout%.txt}_${output_suffix}.txt $contigs "

    # Add minimum-hit-groups if not default
    if [[ "$hit_group" != "default" ]]; then
        kraken_cmd="$kraken_cmd --minimum-hit-groups $hit_group"
    fi

    # Add confidence if not default
    if [[ "$confidence" != "default" ]]; then
        kraken_cmd="$kraken_cmd --confidence $confidence"
    fi

    # Run Kraken2 and extract reads
    echo "Running Kraken2 with --minimum-hit-groups=$hit_group, --confidence=$confidence" >> "$log"
    eval "$kraken_cmd" 2>> "$log"
    
    # Extract classified reads
    awk '$1 == "C" {print $2}' "${krakenout%.txt}_$output_suffix.txt" > "${readslist%.txt}_$output_suffix.txt"
    
    # Count the number of reads
    local num_reads=$(wc -l < "${readslist%.txt}_$output_suffix.txt")

    # Log the result
    echo "Combination: --minimum-hit-groups=$hit_group, --confidence=$confidence, Read Count: $num_reads" >> "$log"

abs_diff_num_reads=$(( num_reads > reads_threshold ? num_reads - reads_threshold : reads_threshold - num_reads ))
abs_diff_closest=$(( closest_read_count > reads_threshold ? closest_read_count - reads_threshold : reads_threshold - closest_read_count ))

# Handle the first iteration or if the current num_reads is closer to the threshold
if (( closest_read_count == 0 || abs_diff_num_reads < abs_diff_closest )); then
    closest_read_count=$num_reads
    closest_conf=$confidence
    closest_hit_groups=$hit_group
    closest_output="${readslist%.txt}_$output_suffix.txt"
    closest_krakenout="${krakenout%.txt}_$output_suffix.txt"
    closest_krakenreport="${krakenreport%.txt}_$output_suffix.txt"
fi
}

# Loop over all combinations of --minimum-hit-groups and --confidence
for hit_group in "${hit_groups[@]}"; do
    for confidence in "${confidences[@]}"; do
        output_suffix="hitgroups_${hit_group}_conf_${confidence}"
        run_kraken "$hit_group" "$confidence" "$output_suffix"
    done
done

# Log the closest combination
echo "Selected combination: --minimum-hit-groups=$closest_hit_groups, --confidence=$closest_conf, Read Count: $closest_read_count" >> "$log"

# If the closest result has more reads than the threshold, randomly subset it
if (( $closest_read_count > $reads_threshold )); then
    echo "Read count exceeds threshold ($reads_threshold). Randomly subsampling to $reads_threshold reads." >> "$log"
    shuf -n "$reads_threshold" "$closest_output" > "$readslist"
else
    # Copy the closest result to the final output file if under or equal to the threshold
    cp "$closest_output" "$readslist"
fi

# Rename the closest Kraken2 output and report to the final names
cp "$closest_krakenout" "$krakenout"
cp "$closest_krakenreport" "$krakenreport"

echo "Final read count in $readslist: $(wc -l < "$readslist")" >> "$log"
