#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_genome_download
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH --time 1:00:00
#SBATCH --partition io


eval "$(conda shell.bash hook)"

conda activate snakemake7

# read in required inputs
while getopts 'i:n:o:p:z:' c
do
  case $c in
    i) input=$OPTARG ;;
    n) API_KEY=$OPTARG ;;
    o) output=$OPTARG ;;
    p) outdir=$OPTARG ;;
    z) programdir=$OPTARG ;;
  esac
done

echo "checking for presence of downloaded genome, removing if exists before redownloading"

rm "${programdir}${outdir}/${output}"

echo "input file : $input"
echo "output file   : $output"
echo "output directory : $outdir"
echo "Program directory   : $programdir"

# api key here, but should make it a variable in the config file and feed it in in future.
export NCBI_API_KEY="$API_KEY"



input2="$programdir""$input"
# Now the actual code to download the genomeS
# 8/4/24 Update to download 1 WGS assembly genome and if present 1 Master TSA transcriptome. 
# (This should hopefully address the huge number of transcript variants making it through from host) 

        sleep $[ ( $RANDOM % 2 )  + 1 ]s

        genus=(`awk '{print $1}' ${input2}`)
        species=(`awk '{print $2}' ${input2}`)
        cd $outdir
        echo " getting genome for: "${genus}" "${species}" "
        echo " Identifying any complete genomes "
        GENOME=(`esearch -db assembly -query ""${genus}" "${species}""[orgn] | efetch -format docsum -stop 2 | tee ""${genus}" "${species}".genome.esearch.docsum"`)
        GENOMEPATHS=(`esearch -db assembly -query ""${genus}" "${species}""[orgn] | efetch -format docsum -stop 2 | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee ""${genus}" "${species}".genome.esearch.docsum"`)

if [ ! ${#GENOMEPATHS[@]} -eq 0 ]; then



for element in "${GENOMEPATHS[@]}"; do
    # Check if the element is not already in the unique array
    if [[ ! " ${GENOMEPATHSuniq[@]} " =~ " ${element} " ]]; then
        # Add the unique element to the unique array
        GENOMEPATHSuniq+=("$element")
    fi
done

for element in "${GENOMEPATHSuniq[@]}"; do
    # Check if the element does not start with a space
    if [[ "$element" != " "* ]]; then
        # Add the element to the new array
        GENOMEPATHSuniq2+=("$element")
    fi
done




START=0
END=$((${#GENOMEPATHSuniq2[@]} - 1))

for (( c=$START; c<=$END; c++ ))
do

echo "iteration $c"

GENOMEPATH2=("${GENOMEPATHSuniq2[c]}")

GENOMEPATH2=${GENOMEPATH2//ftp:/https:}

genomebase=(`basename $GENOMEPATH2`)

echo "${genomebase[0]}"
echo "${GENOMEPATH2[0]}"


sleep $[ ( $RANDOM % 1 )  + 1 ]s

echo " Now getting genome "${programdir}${outdir}/${output}.${c}" which is found at "${GENOMEPATH2[0]}"/"${genomebase[0]}"_genomic.fna.gz "

wget -O "${programdir}${outdir}/${output}.${c}" "${GENOMEPATH2[0]}"/"${genomebase[0]}"_genomic.fna.gz

zcat "${programdir}${outdir}/${output}.${c}" >> "${programdir}${outdir}/${output}"
rm "${programdir}${outdir}/${output}.${c}"

done

seqkit rmdup -n "${programdir}${outdir}/${output}_predup" -o "${programdir}${outdir}/${output}"


rm "${programdir}${outdir}/${output}_predup"


fi


if [ ${#GENOMEPATHS[@]} -eq 0 ]; then

echo " no genomes were found for the target species. Moving to the target genus. Up to four genomes will be downloaded randomly from genomes within the identified genus "



GENOME=(`esearch -db assembly -query "${genus}"[orgn] | efetch -format docsum -stop 4 | tee ""${genus}" "${species}".genome.esearch.docsum"`)
GENOMEPATHS=(`esearch -db assembly -query "${genus}"[orgn] | efetch -format docsum -stop 4 | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee ""${genus}" "${species}".genome.esearch.docsum"`)


if [ ! ${#GENOMEPATHS[@]} -eq 0 ]; then

for element in "${GENOMEPATHS[@]}"; do
    # Check if the element is not already in the unique array
    if [[ ! " ${GENOMEPATHSuniq[@]} " =~ " ${element} " ]]; then
        # Add the unique element to the unique array
        GENOMEPATHSuniq+=("$element")
    fi
done

for element in "${GENOMEPATHSuniq[@]}"; do
    # Check if the element does not start with a space
    if [[ "$element" != " "* ]]; then
        # Add the element to the new array
        GENOMEPATHSuniq2+=("$element")
    fi
done


START=0
END=$((${#GENOMEPATHSuniq2[@]} - 1))






START=0
END=$((${#GENOMEPATHS[@]} - 1))

for (( c=$START; c<=$END; c++ ))
do

echo "iteration $c"

GENOMEPATH2=("${GENOMEPATHSuniq2[c]}")

GENOMEPATH2=${GENOMEPATH2//ftp:/https:}

genomebase=(`basename $GENOMEPATH2`)

echo "${genomebase[0]}"
echo "${GENOMEPATH2[0]}"


sleep $[ ( $RANDOM % 1 )  + 1 ]s

echo " getting "${programdir}${outdir}/${output}.${c}" "

wget -O "${programdir}${outdir}/${output}.${c}" "${GENOMEPATH2[0]}"/"${genomebase[0]}"_genomic.fna.gz

zcat "${programdir}${outdir}/${output}.${c}" >> "${programdir}${outdir}/${output}_predup"
rm "${programdir}${outdir}/${output}.${c}"

done
seqkit rmdup -n "${programdir}${outdir}/${output}_predup" -o "${programdir}${outdir}/${output}"


rm "${programdir}${outdir}/${output}_predup"
fi

fi

# Figure out where this needs to save for while loop to be added! 

touch finished.txt

exit