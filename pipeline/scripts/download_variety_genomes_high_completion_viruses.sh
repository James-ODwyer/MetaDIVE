#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_genome_download_multi
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 4
#SBATCH --mem 12G
#SBATCH --time 24:00:00
#SBATCH --partition io


eval "$(conda shell.bash hook)"

conda activate snakemake7

while getopts 'g:j:k:n:p:r:t:z:' c
do
  case $c in
    g) inputgenomedir=$OPTARG ;;
    j) RawSSUreads1=$OPTARG ;;
    k) RawSSUreads2=$OPTARG ;;
    n) API_KEY=$OPTARG ;;
    p) inputhostremovedir=$OPTARG ;;
    r) outputconsensusdir=$OPTARG ;;
    t) samplename=$OPTARG ;;
    z) programdir=$OPTARG ;;
  esac
done

echo " input genome dir $inputgenomedir "
echo " raw CO1 reads 1 $RawSSUreads1 "
echo " raw CO1 reads 2 $RawSSUreads2 "
echo " input host removed dir $inputhostremovedir "
echo " output consensus dir $outputconsensusdir "
echo " sample name: $samplename "
echo " program dir: $programdir "


export NCBI_API_KEY='$API_KEY'


# Edited version. New pathway takes taxid's from identified genomes with minimum 100 hits to the genome. 


# list taxids to run

cd "${programdir}${inputgenomedir}"

numgenomes=(`ls *TaxId* | wc -l`)

if [ ${numgenomes} -ge 1 ]; then


genometaxids=(`ls *TaxId*`)

for file in ${genometaxids[@]}; do
cd "${programdir}${inputgenomedir}"

echo "file $file"

speciesname=(`awk '{sub(/_TaxId.*/,""); print}' <<<$file`)

# take top taxid code from file (all will be identical within each document)
read -r taxid < "$file"

sleep $[ ( $RANDOM % 2 )  + 5 ]s
cd "${programdir}${outputconsensusdir}"

# First part. Generate a docsum for sequence length for assemblies. 
# Will return the length of all assemblies. 

genomelengths=(`esearch -db assembly -query "${taxid}"[TaxId] | elink -target nuccore | efetch -format docsum | xtract -pattern DocumentSummary -element Slen`)

sleep $[ ( $RANDOM % 2 )  + 5 ]s

Id=(`esearch -db taxonomy -query "${taxid}[TaxId]" | elink -target nuccore | efetch -format docsum -start 1 -stop 1 | xtract -pattern DocumentSummary -element Organism`)
orgname1=(`echo "${Id[@]}"`)
orgname2="${orgname1[@]}"
orgname=(" ${orgname2} ")

if [[ "${orgname[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname[@]:1}")
fi
if [[ "${orgname[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname[@]:1}")
fi

if [ ${#genomelengths[@]} -gt 0 ]; then

echo " At least one whole genome was found in the assembly database, now generating average genome length \n and downloading genomes at minimum 50% size "

avggenomeength=(`echo "${genomelengths[@]}" | tr ' ' '\n' | awk '{sum += $1} END {print sum / NR}'`)

echo " average genome length for geones of $file $avggenomeength bases "

min_return_length=$(echo "$avggenomeength / 4" | bc -l)
max_return_length=$(echo "$avggenomeength * 4" | bc -l)

sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db nuccore -query "${orgname}[Organism] AND (${min_return_length}[SLEN] : ${max_return_length}[SLEN])" | efetch -format fasta -stop 50 > ${speciesname}_combined_references.fasta

full_genomes=(`grep ">" ${speciesname}_combined_references.fasta | wc -l`)
if [ ${full_genomes} -le 15 ]; then 

echo " fewer than 15 genomes found, expanding size cut off to 12% of average genome length to allow fragments "
min_return_length=$(echo "$avggenomeength / 8" | bc -l)

sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db nuccore -query "${orgname}[Organism] AND ${min_return_length}:${max_return_length}[SLEN]" | efetch -format fasta -stop 50 > ${speciesname}_combined_references.fasta

fi


fi

if [ ${#genomelengths[@]} -eq 0 ]; then

echo " No assembled viral genomes were found, now searching nuccore and will assume the largest avaliable sequence is near completion"
sleep $[ ( $RANDOM % 2 )  + 5 ]s
genomelengths=(`esearch -db nuccore -query "${orgname}[Organism]" | efetch -format docsum | xtract -pattern DocumentSummary -element Slen`)


max=${genomelengths[0]}

for value in "${genomelengths[@]}"
do
    if (( value > max )); then
        max=$value
    fi
done

min_return_length=$(echo "$max / 8" | bc -l)

sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db nuccore -query "${orgname}[Organism] AND ${min_return_length}:9999999[SLEN]" | efetch -format fasta -stop 50 > ${speciesname}_combined_references.fasta



full_genomes=(`grep ">" ${speciesname}_combined_references.fasta | wc -l`)


if [ ${full_genomes} -le 15 ]; then 

echo " Fewer than 15 sequences were found, changing minimum size from 12% of max length sequence to 5% "
min_return_length=$(echo "$max/ 20" | bc -l)

sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db nuccore -query "${orgname}[Organism] AND ${min_return_length}[SLEN]:9999999" | efetch -format fasta -stop 50 > ${speciesname}_combined_references.fasta
fi



fi

sleep $[ ( $RANDOM % 2 )  + 1 ]s


workdir=`pwd`

if [ ! -d "${workdir}/reference_genomes" ]; then

mkdir "${workdir}/reference_genomes"
fi

R1=${programdir}${RawSSUreads1}
R2=${programdir}${RawSSUreads2}

if [ ! -d "${workdir}/reference_genomes/${speciesname}" ]; then
mkdir "${workdir}/reference_genomes/${speciesname}"
fi

if [ ! -d "${workdir}/hits_${speciesname}" ]; then
mkdir "${workdir}/hits_${speciesname}"
fi

echo " running bowtie on genomes "
bowtie2-build ${speciesname}_combined_references.fasta "${workdir}/reference_genomes/${speciesname}/index_genome" --threads ${SLURM_CPUS_PER_TASK} \
            --large-index && \
        bowtie2 -x "${workdir}/reference_genomes/${speciesname}/index_genome" -1 ${R1} -2 ${R2} \
            -p ${SLURM_CPUS_PER_TASK} \
            -S "${workdir}/hits_${speciesname}/${speciesname}_sam.sam" \
            --fast


samtools view -@ ${SLURM_CPUS_PER_TASK} -F 4 -h ${workdir}/hits_${speciesname}/${speciesname}_sam.sam > ${workdir}/hits_${speciesname}/${speciesname}_sam_hits.sam

samtools fastq -@ ${SLURM_CPUS_PER_TASK} --verbosity 1 -n ${workdir}/hits_${speciesname}/${speciesname}_sam_hits.sam \
-1 "${workdir}/hits_${speciesname}/${speciesname}_hits_R1.fastq" -2 "${workdir}/hits_${speciesname}/${speciesname}_hits_R2.fastq"


line_count=$(wc -l < "${workdir}/hits_${speciesname}/${speciesname}_hits_R1.fastq")

# Check if the line count is greater than 2
if [ "$line_count" -gt 100 ]; then

echo "multiple reads were aligned using bowtie. Now generating spades contigs"

samtools view ${workdir}/hits_${speciesname}/${speciesname}_sam_hits.sam | awk '{counts[$3]++} END {for (genome in counts) print genome, counts[genome]}' > ${speciesname}_reads_per_genome.txt
genomehits_per_ref=(`samtools view ${workdir}/hits_${speciesname}/${speciesname}_sam_hits.sam | awk '{counts[$3]++} END {for (genome in counts) print genome, counts[genome]}'`)


topaccessionhit=(`awk 'BEGIN {max_value = 0; max_name = ""}
     $2 > max_value {max_value = $2; max_name = $1}
     END {print max_name}' ${speciesname}_reads_per_genome.txt`)

echo " the closest matching genome to the returned virus in this study was $topaccessionhit "

awk -v target="${topaccessionhit}" '/^>/{if (name ~ target) print name ORS seq; name=$0; seq=""; next} {seq = seq $0} END {if (name ~ target) print name ORS seq}' ${speciesname}_combined_references.fasta > ${topaccessionhit}_${speciesname}.fasta


if [ ! -d "${workdir}/reference_genomes/${speciesname}/closest_relative" ]; then
mkdir "${workdir}/reference_genomes/${speciesname}/closest_relative"
fi

if [ ! -d "${workdir}/hits_${speciesname}/closest_relative" ]; then
mkdir "${workdir}/hits_${speciesname}/closest_relative"
fi



bowtie2-build ${topaccessionhit}_${speciesname}.fasta "${workdir}/reference_genomes/${speciesname}/closest_relative/index_genome" --threads ${SLURM_CPUS_PER_TASK} \
            --large-index && \
        bowtie2 -x "${workdir}/reference_genomes/${speciesname}/closest_relative/index_genome" -1 ${R1} -2 ${R2} \
            --al-conc "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_R%.fastq" \
            -p ${SLURM_CPUS_PER_TASK} \
            -S "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam" \
            --fast

samtools view -@ ${SLURM_CPUS_PER_TASK} -F 4 -h ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam > ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam_hits.sam

samtools fastq -@ ${SLURM_CPUS_PER_TASK} --verbosity 1 -n ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam_hits.sam \
-1 "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_R1.fastq" -2 "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_R2.fastq" -U "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_singletons.fastq"




samtools view -S -b -F 4 -h ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam > ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits.bam

samtools sort ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits.bam -o ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_sorted.bam

        samtools consensus -d 3 -a -m "simple" -c 0.51 --use-qual --show-ins yes --show-del no ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_hits_sorted.bam -o ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta



echo " running Spades on all bowtie hits "


conda activate spades

spades.py \
            -1 "${workdir}/hits_${speciesname}/${speciesname}_hits_R1.fastq" \
            -2 "${workdir}/hits_${speciesname}/${speciesname}_hits_R2.fastq" \
            --isolate \
            -t ${SLURM_CPUS_PER_TASK} \
            -m 12 \
            -o "${workdir}/hits_${speciesname}/spades"


conda deactivate

Missingdata=`grep --no-group-separator -o -i 'N' ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta | wc -l`


fastasize=`wc -m ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta`



echo " The closest relative alignment for ${speciesname} shows ${Missingdata} Ns in an estimated genome size of ${fastasize} bases"

proportion_N=$(echo "scale=4; ${Missingdata} / $(echo ${fastasize} | awk '{print $1}')" | bc)


if (( $(echo "${proportion_N} < 0.2500" | bc -l) )); then



if [ ! -d "${workdir}/finished_genomes" ]; then
mkdir "${workdir}/finished_genomes"
fi

if [ ! -d "${workdir}/finished_genomes/${speciesname}" ]; then
mkdir "${workdir}/finished_genomes/${speciesname}"
fi




cp ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta ${workdir}/finished_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta
cp ${workdir}/hits_${speciesname}/spades/scaffolds.fasta ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta


for file in ${workdir}/hits_${speciesname}/*R1*; do
  cp "$file" "${workdir}/finished_genomes/"
  gzip "${workdir}/finished_genomes/$(basename "$file")"
done

for file in ${workdir}/hits_${speciesname}/*R2*; do
  cp "$file" "${workdir}/finished_genomes/"
  gzip "${workdir}/finished_genomes/$(basename "$file")"
done


sed -i 's/^>.*$/>'"${samplename}_aligned_to_${topaccessionhit}"'/' ${workdir}/finished_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta

rm "${workdir}/hits_${speciesname}/${speciesname}_sam.sam"
rm "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam"

NSEQS=2

awk "/^>/ {n++} n>$NSEQS {exit} {print}" ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta > ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_1_length.*$/>'"${samplename}_scaffold1"'/' ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_2_length.*$/>'"${samplename}_scaffold2"'/' ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta

if [ ! -d "${workdir}/finished_genomes/${speciesname}_multifasta" ]; then
mkdir "${workdir}/finished_genomes/${speciesname}_multifasta"
fi

cat ${workdir}/finished_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta ${workdir}/finished_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta ${speciesname}_combined_references.fasta > ${workdir}/finished_genomes/${speciesname}_multifasta/${speciesname}_combined_fasta.fasta

fi

# repeat for if sequence is middling missing data
if (( $(echo "${proportion_N} > 0.2500 && ${proportion_N} < 0.7500" | bc -l) )); then

if [ ! -d "${workdir}/moderate_coverage_genomes" ]; then
mkdir "${workdir}/moderate_coverage_genomes"
fi

if [ ! -d "${workdir}/moderate_coverage_genomes/${speciesname}" ]; then
mkdir "${workdir}/moderate_coverage_genomes/${speciesname}"
fi

cp ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta
cp ${workdir}/hits_${speciesname}/spades/scaffolds.fasta ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta

for file in ${workdir}/hits_${speciesname}/*R1*; do
  cp "$file" "${workdir}/moderate_coverage_genomes/"
  gzip "${workdir}/moderate_coverage_genomes/$(basename "$file")"
done

for file in ${workdir}/hits_${speciesname}/*R2*; do
  cp "$file" "${workdir}/moderate_coverage_genomes/"
  gzip "${workdir}/moderate_coverage_genomes/$(basename "$file")"
done

sed -i 's/^>.*$/>'"${samplename}_aligned_to_${topaccessionhit}"'/' ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta

rm "${workdir}/hits_${speciesname}/${speciesname}_sam.sam"
rm "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam"

NSEQS=2

awk "/^>/ {n++} n>$NSEQS {exit} {print}" ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta > ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_1_length.*$/>'"${samplename}_scaffold1"'/' ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_2_length.*$/>'"${samplename}_scaffold2"'/' ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta

if [ ! -d "${workdir}/moderate_coverage_genomes/${speciesname}_multifasta" ]; then
mkdir "${workdir}/moderate_coverage_genomes/${speciesname}_multifasta"
fi

cat ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta ${workdir}/moderate_coverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta ${speciesname}_combined_references.fasta > ${workdir}/moderate_coverage_genomes/${speciesname}_multifasta/${speciesname}_combined_fasta.fasta

fi








# repeat for if the majority of sequence is N 
if (( $(echo "${proportion_N} > 0.7500" | bc -l) )); then


if [ ! -d "${workdir}/lowcoverage_genomes" ]; then
mkdir "${workdir}/lowcoverage_genomes"
fi

if [ ! -d "${workdir}/lowcoverage_genomes/${speciesname}" ]; then
mkdir "${workdir}/lowcoverage_genomes/${speciesname}"
fi

cp ${workdir}/hits_${speciesname}/closest_relative/${speciesname}_consensus_closest_relative.fasta ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta
cp ${workdir}/hits_${speciesname}/spades/scaffolds.fasta ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta

sed -i 's/^>.*$/>'"${samplename}_aligned_to_${topaccessionhit}"'/' ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta

rm "${workdir}/hits_${speciesname}/${speciesname}_sam.sam"
rm "${workdir}/hits_${speciesname}/closest_relative/${speciesname}_sam.sam"

NSEQS=2

awk "/^>/ {n++} n>$NSEQS {exit} {print}" ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo.fasta > ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_1_length.*$/>'"${samplename}_scaffold1"'/' ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta
sed -i 's/^>NODE_2_length.*$/>'"${samplename}_scaffold2"'/' ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta

if [ ! -d "${workdir}/lowcoverage_genomes/${speciesname}_multifasta" ]; then
mkdir "${workdir}/lowcoverage_genomes/${speciesname}_multifasta"
fi

cat ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_generated_scaffolds_denovo_top2.fasta ${workdir}/lowcoverage_genomes/${speciesname}/${speciesname}_consensus_sequence_from_closest_relative.fasta ${speciesname}_combined_references.fasta > ${workdir}/lowcoverage_genomes/${speciesname}_multifasta/${speciesname}_combined_fasta.fasta

fi

fi 


if [ "$line_count" -le 20 ]; then

touch "${programdir}${outputconsensusdir}/${speciesname}_failed_minimum_read_threshold"





fi

done

# This ends with the following results per genome, 
# 1. fasta of all good genomes for references
# 2. spades assembly per genome for reads within this study.
# 3. a combined fasta ready for alignment

# need new script to blast top scaffold for consensus mapping
# them align and then finally trimming and tree


touch "${programdir}${outputconsensusdir}/finished2.txt"


fi


if [ ${numgenomes} -lt 1 ]; then

touch "${programdir}${outputconsensusdir}/finished2.txt"

fi
