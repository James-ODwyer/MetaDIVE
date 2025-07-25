#!/bin/bash

set -euo pipefail

# === DEFAULTS ===
SCRIPT_LOCATIONS=("envs" "pipeline/scripts" "pipeline" "databases")
HEADER_PREFIX=""
BATCH_STYLE=""
PARTITION=""
ACCOUNT=""
DOWNLOAD_PARTITION="none"
NODES_SUPPORTED="yes"
NTASKS_SUPPORTED="yes"
CPUS_SUPPORTED="yes"

print_usage() {
    echo "Usage: $0 --batch-system SYSTEM --partition NAME [--account NAME] [--download_partition NAME]"
    echo "             [--nodes-supported yes|no] [--ntasks-per-node-supported yes|no] [--cpus-per-task-supported yes|no]"
    echo ""
    echo "Supported batch systems: SLURM, PBS, QSUB, LSF, SGE, CUSTOM"
}

# === ARG PARSER ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-system)
            BATCH_SYSTEM=$(echo "$2" | tr '[:lower:]' '[:upper:]')
            case $BATCH_SYSTEM in
                SLURM) HEADER_PREFIX="#SBATCH"; BATCH_STYLE="slurm";;
                PBS|QSUB) HEADER_PREFIX="#PBS"; BATCH_STYLE="pbs";;
                LSF) HEADER_PREFIX="#BSUB"; BATCH_STYLE="lsf";;
                SGE) HEADER_PREFIX="#$"; BATCH_STYLE="sge";;
                CUSTOM) HEADER_PREFIX="#CUSTOM"; BATCH_STYLE="custom";;
                *) echo "Unsupported batch system: $2"; print_usage; exit 1;;
            esac
            shift 2;;
        --partition) PARTITION="$2"; shift 2;;
        --account) ACCOUNT="$2"; shift 2;;
        --download_partition) DOWNLOAD_PARTITION="$2"; shift 2;;
        --nodes-supported) NODES_SUPPORTED="$2"; shift 2;;
        --ntasks-per-node-supported) NTASKS_SUPPORTED="$2"; shift 2;;
        --cpus-per-task-supported) CPUS_SUPPORTED="$2"; shift 2;;
        -h|--help) print_usage; exit 0;;
        *) echo "Unknown option: $1"; print_usage; exit 1;;
    esac
done

if [[ -z "$BATCH_STYLE" || -z "$PARTITION" ]]; then
    echo "Error: Missing required arguments."
    print_usage
    exit 1
fi

# === PARAMETER NAME MAPPING ===
replace_param_names() {
    local line="$1"
    local is_download_script="$2"
    local new_partition

    # Strip any known prefix (e.g., #SBATCH, #PBS)
    line=$(echo "$line" | sed -E 's/^#(SBATCH|PBS|BSUB|\$|CUSTOM)[[:space:]]+//')

    # Check if line should be removed due to unsupported resource flag
    if [[ "$line" =~ --nodes[[:space:]=] ]] && [[ "$NODES_SUPPORTED" == "no" ]]; then
        return 1
    fi
    if [[ "$line" =~ --ntasks-per-node[[:space:]=] ]] && [[ "$NTASKS_SUPPORTED" == "no" ]]; then
        return 1
    fi
    if [[ "$line" =~ --cpus-per-task[[:space:]=] ]] && [[ "$CPUS_SUPPORTED" == "no" ]]; then
        return 1
    fi

    # Decide which partition to use
    if [[ "$DOWNLOAD_PARTITION" != "none" && "$is_download_script" == "true" ]]; then
        new_partition="$DOWNLOAD_PARTITION"
    else
        new_partition="$PARTITION"
    fi

    # Rewrite resource parameter names for each batch system
    case $BATCH_STYLE in
        slurm)
            line=$(echo "$line" | sed -E "s/--partition[ =][^ ]*/--partition $new_partition/")
            [[ -n $ACCOUNT ]] && line=$(echo "$line" | sed -E "s/--account[ =][^ ]*/--account=$ACCOUNT/")
            ;;
        pbs)
            line=$(echo "$line" | sed -E "s/--cpus-per-task[ =]([0-9]+)/-l nodes=1:ppn=\1/")
            line=$(echo "$line" | sed -E "s/--mem[ =]([0-9A-Za-z]+)/-l mem=\1/")
            line=$(echo "$line" | sed -E "s/--partition[ =][^ ]*/-q $new_partition/")
            [[ -n $ACCOUNT ]] && line=$(echo "$line" | sed -E "s/--account[ =][^ ]*/-A $ACCOUNT/")
            ;;
        lsf)
            line=$(echo "$line" | sed -E "s/--cpus-per-task[ =]([0-9]+)/-n \1/")
            line=$(echo "$line" | sed -E "s/--mem[ =]([0-9A-Za-z]+)/-R \"rusage[mem=\1]\"/")
            line=$(echo "$line" | sed -E "s/--partition[ =][^ ]*/-q $new_partition/")
            [[ -n $ACCOUNT ]] && line=$(echo "$line" | sed -E "s/--account[ =][^ ]*/-P $ACCOUNT/")
            ;;
        sge)
            line=$(echo "$line" | sed -E "s/--cpus-per-task[ =]([0-9]+)/-pe smp \1/")
            line=$(echo "$line" | sed -E "s/--mem[ =]([0-9A-Za-z]+)/-l mem_free=\1/")
            line=$(echo "$line" | sed -E "s/--partition[ =][^ ]*/-q $new_partition/")
            [[ -n $ACCOUNT ]] && line=$(echo "$line" | sed -E "s/--account[ =][^ ]*/-A $ACCOUNT/")
            ;;
        custom)
            line="$line"  # No changes for custom
            ;;
    esac

    # Final check: conditionally remove partition line
    if [[ "$line" =~ --partition[[:space:]=] || "$line" =~ -q[[:space:]] ]]; then
        if [[ "$is_download_script" == "true" && "$DOWNLOAD_PARTITION" != "none" ]]; then
            :
        elif [[ "$PARTITION" == "none" || -z "$PARTITION" ]]; then
            return 1
        fi
    fi

    echo "$line"
}

# === MODIFY FILES ===
for DIR in "${SCRIPT_LOCATIONS[@]}"; do
    [[ ! -d "$DIR" ]] && continue
    find "$DIR" -type f -name "*.sh" | while read -r FILE; do
        echo "Modifying $FILE..."
        cp "$FILE" "$FILE.bak"

        head -n 10 "$FILE" > tmp_head.txt
        tail -n +11 "$FILE" > tmp_tail.txt

        if [[ "$DOWNLOAD_PARTITION" != "none" ]]; then
            if grep -q -- "--partition io" tmp_head.txt; then
                script_contains_download_io="true"
            else
                script_contains_download_io="false"
            fi
        else
            script_contains_download_io="false"
        fi

        > tmp_head_fixed.txt
        while read -r line; do
            if [[ "$line" =~ ^#(SBATCH|PBS|BSUB|\$|CUSTOM)[[:space:]]+ ]]; then
                if fixed_line=$(replace_param_names "$line" "$script_contains_download_io"); then
                    echo "$HEADER_PREFIX $fixed_line" >> tmp_head_fixed.txt
                else
                    echo "Removed unsupported line in $FILE: $line"
                fi
            else
                echo "$line" >> tmp_head_fixed.txt
            fi
        done < tmp_head.txt

        cat tmp_head_fixed.txt tmp_tail.txt > "$FILE"
        rm tmp_head.txt tmp_tail.txt tmp_head_fixed.txt
    done
done

echo -e "\n Cleaning up backup files..."
find "${SCRIPT_LOCATIONS[@]}" -type f -name "*.sh.bak" -exec rm -f {} \;

echo -e "\n Batch headers updated for batch system: $BATCH_SYSTEM"
