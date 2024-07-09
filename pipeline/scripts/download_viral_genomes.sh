#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_genome_download
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 2:00:00
#SBATCH --partition io

eval "$(conda shell.bash hook)"

conda activate snakemake7


# read in required inputs
# input is the annotated gene file from Dramv and Matt's py code to combine the annot file with the gff gene file
# outdir will be the 16_REFERENCE_VIRAL_GENOMES/{sample}/
# programdir will be the config["programdir"]


while getopts 'j:n:p:z:' c
do
  case $c in
    j) input2=$OPTARG ;;
    n) API_KEY=$OPTARG ;;
    p) outdir=$OPTARG ;;
    z) programdir=$OPTARG ;;
  esac
done


export NCBI_API_KEY='$API_KEY'


# Input 2 is the taxid file. Need to update rule still with the argpass call.

echo "input file of annotated genes: $input2"
echo "output directory : $outdir"
echo "Program directory   : $programdir"




declare -a species=()
declare -a genus=()
declare -a orgname=()
declare -a orgtax=()
# now to add input 2 into the lists 

taxid_values=()

# Read the TSV file line by line and extract the first column value
while IFS=$'\t' read -r first_col _; do
    taxid_values+=("$first_col")
done < "$input2"

taxid_values=("${taxid_values[@]% }")

for id in ${taxid_values[@]}; do

sleep $[ ( $RANDOM % 5 ) ]s

value="$id"
value="${value#"${value%%[![:space:]]*}"}"
value="${value%"${value##*[![:space:]]}"}"

echo "$value"


#break retry loop in case NCBI is busy (for Organism)
retry_count=0
max_retries=5
retry_wait_time=120  # 120 seconds
while true; do
    # Run the eutilities command
    Id=$(esearch -db taxonomy -query "$value"[TaxId] | elink -target nuccore | efetch -format docsum -start 1 -stop 1 | xtract -pattern DocumentSummary -element Organism)

    # Check if the command was successful
 if [[ -n "$Id" && ! "$Id" =~ ^[[:space:]]*$ ]]; then
        echo "Command executed successfully"
        break
    else
        echo "Error executing the organism collection command. Retrying..."

        # Increment retry count
        ((retry_count++))
        echo " retry count at $retry_count "

        # Check if max retries exceeded
        if [ $retry_count -ge $max_retries ]; then
            echo "Max retries exceeded. OrgID is assigned as Unknown to preserve order of sequences. Exiting..."
            Id="Unknown"
            echo "species name $Id "
            break
        fi

        # Wait for some time before retrying
        echo "Waiting for $retry_wait_time seconds before retrying..."
        sleep $retry_wait_time
    fi
done

sleep $[ ( $RANDOM % 1 ) + 5 ]s
    # Run the eutilities command
    Idtax="$value"

#echo "organism name ${Idnodup[@]}"

genus1=(`echo "${Id[0]}"`)
species1=(`echo "${Id[@]}"`)
orgname1=(`echo "${Id[@]}"`)
orgname2="${orgname1[@]}"


if [[ "${Idtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${Idtax[@]:1}")
fi
if [[ "${Idtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${Idtax[@]:1}")
fi
if [[ "${Idtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${Idtax[@]:1}")
fi
if [[ "${orgname2[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname2[@]:1}")
fi
if [[ "${orgname2[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname2[@]:1}")
fi
if [[ "${orgname2[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgname2[@]:1}")
fi

if ! grep -qw "${Idtax[0]}" <<< "${orgtax[*]}"; then
echo " Novel species to be added" 

echo "${Idtax[0]}" 
echo "${species1[@]}" 
genus+=(${genus1})
species+=("${species1[@]}")
orgname+=(" ${orgname2} ")
orgtax+=("${Idtax}")


orgname2=""
species1=""
genus1=""

fi

echo "finished $value "

done


START=0
#END=(${#orgname[@]})
END=($((${#orgname[@]}-1)))

echo $START
echo $END

sleep $[ ( $RANDOM % 2 )  + 1 ]s

cd $outdir

pwd

# remove excess spaces potentially front loaded onto taxids codes
if [[ "${orgname[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname[@]:1}")
fi
if [[ "${orgname[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgname=("${orgname[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi
if [[ "${orgtax[0]}" =~ ^[[:space:]]*$ ]]; then
  # Remove the first element using array slicing
  orgtax=("${orgtax[@]:1}")
fi


for (( c=$START; c<=$END; c++ ))
do

orgnameid=("${orgname[$c]}")
orgnameidnospace=$(echo -n "${orgnameid//[[:space:].,:;-&]/}")


orgtaxsingle=("${orgtax[$c]}")

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"

echo " Running organism $orgnameidnospace with taxID of $orgtaxsingle "

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 1 ]s

echo " Running esearch to generate docusum for ${orgnameid} "


echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 3 ]s


GENOME=`esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1  | xtract -pattern DocumentSummary -element FtpPath_RefSeq | tee "${orgnameidnospace}.genome.esearch.docsum"`

sleep $[ ( $RANDOM % 2 )  + 5 ]s

esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | tee "${orgnameidnospace}_wholedocsumdata.docsum"
sleep $[ ( $RANDOM % 2 )  + 5 ]s


if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ]
then

GENOME=`esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1  | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee "${orgnameidnospace}.genome.esearch.docsum"`

sleep $[ ( $RANDOM % 2 )  + 5 ]s

esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | tee "${orgnameidnospace}_wholedocsumdata.docsum"
sleep $[ ( $RANDOM % 2 )  + 5 ]s

fi


if [ -s "${orgnameidnospace}.genome.esearch.docsum" ] 
then
echo "${GENOME[@]}"

GENOME2=("${GENOME[0]}")

GENOME2=${GENOME2//ftp:/https:}

genomebase=(`basename $GENOME2`)
echo "${genomebase[0]}"
echo "${GENOME2[0]}"

sleep $[ ( $RANDOM % 2 )  + 5 ]s

wget -O "${orgnameidnospace}_genome.fna.gz" "${GENOME2[0]}"/"${genomebase[0]}"_genomic.fna.gz

sleep $[ ( $RANDOM % 2 )  + 5 ]s

echo " genome document summary ${GENOME2[0]}     ( ${orgnameid} ) "

fi


if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ]; then
sleep $[ ( $RANDOM % 2 )  + 5 ]s

echo "No RefSeq genome identified, now trying genome database "

GENOME=`esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee "${orgnameidnospace}.genome.esearch.docsum"`
sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"

echo "Running esearch to save top filepath "

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




if [ -s "${orgnameidnospace}.genome.esearch.docsum" ];then


echo "${GENOME[@]}"

GENOME2=("${GENOME[0]}")

GENOME2=${GENOME2//ftp:/https:}

genomebase=(`basename $GENOME2`)
echo "${genomebase[0]}"
echo "${GENOME2[0]}"

wget -O "${orgnameidnospace}_genome.fna.gz" "${GENOME2[0]}"/"${genomebase[0]}"_genomic.fna.gz

echo " genome document summary ${GENOME2[0]}     ( ${orgnameid} ) "


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




fi
fi


if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ]; then
sleep $[ ( $RANDOM % 2 )  + 5 ]s

echo "No genome identified in genome database, now trying assembly database for non refseq "

GENOME=`esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 2 | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee "${orgnameidnospace}.genome.esearch.docsum"`
sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"

echo "Running esearch to save top filepath "

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




if [ -s "${orgnameidnospace}.genome.esearch.docsum" ];then


echo "${GENOME[@]}"

GENOME2=("${GENOME[0]}")

GENOME2=${GENOME2//ftp:/https:}

genomebase=(`basename $GENOME2`)
echo "${genomebase[0]}"
echo "${GENOME2[0]}"

wget -O "${orgnameidnospace}_genome.fna.gz" "${GENOME2[0]}"/"${genomebase[0]}"_genomic.fna.gz

echo " genome document summary ${GENOME2[0]}     ( ${orgnameid} ) "


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"


fi
fi



if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 5 ]s

echo "No complete genomes identified in either the genomes or assembly databases, now trying any genome for                  ${orgnameid}"

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR genome[All Fields]) AND 2000:9999999[Slen]" | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"
echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 5 ]s

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR genome[All Fields]) AND 2000:9999999[Slen]" | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 5 ]s



echo " genome document summary  $GENOME                  ( ${orgnameid} ) "

echo " ${GENOME[@]} "

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




fi

if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 5 ]s

echo "No RefSeq genome identified, now trying any genome                  ( ${orgnameid} ) "

sleep $[ ( $RANDOM % 2 )  + 5 ]s
esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR partial[All Fields]) AND 2000:9999999[Slen]" | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 5 ]s

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR partial[All Fields]) AND 2000:9999999[Slen]" | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"

echo " genome document summary $GENOME                  ( ${orgnameid} ) "

echo " ${GENOME[@]} "


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"


fi

# Redo with wait time in case it errored 
if [ ! -s "${orgnameidnospace}_genome.fna.gz" ];then


echo " Running organism $orgnameidnospace with taxID of $orgtaxsingle . This is a repeat genome download attempt after initial failure. "

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 120 ]s

echo " Running esearch to generate docusum for ${orgnameid} "


echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 3 ]s


GENOME=`esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1  | xtract -pattern DocumentSummary -element FtpPath_RefSeq | tee "${orgnameidnospace}.genome.esearch.docsum"`

sleep $[ ( $RANDOM % 2 )  + 3 ]s

esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | tee "${orgnameidnospace}_wholedocsumdata.docsum"
sleep $[ ( $RANDOM % 2 )  + 3 ]s


if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ]
then

GENOME=`esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1  | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee "${orgnameidnospace}.genome.esearch.docsum"`

sleep $[ ( $RANDOM % 2 )  + 3 ]s

esearch -db assembly -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | tee "${orgnameidnospace}_wholedocsumdata.docsum"
sleep $[ ( $RANDOM % 2 )  + 3 ]s

fi


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"


if [ -s "${orgnameidnospace}.genome.esearch.docsum" ] 
then
echo "${GENOME[@]}"

GENOME2=("${GENOME[0]}")

GENOME2=${GENOME2//ftp:/https:}

genomebase=(`basename $GENOME2`)
echo "${genomebase[0]}"
echo "${GENOME2[0]}"

sleep $[ ( $RANDOM % 2 )  + 3 ]s

wget -O "${orgnameidnospace}_genome.fna.gz" "${GENOME2[0]}"/"${genomebase[0]}"_genomic.fna.gz

sleep $[ ( $RANDOM % 2 )  + 3 ]s

echo " genome document summary ${GENOME2[0]}     ( ${orgnameid} ) "

fi
fi

if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ]; then
sleep $[ ( $RANDOM % 2 )  + 3 ]s

echo "No RefSeq genome identified, now trying genome database "

GENOME=`esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum -stop 1 | xtract -pattern DocumentSummary -element FtpPath_GenBank | tee "${orgnameidnospace}.genome.esearch.docsum"`
sleep $[ ( $RANDOM % 2 )  + 3 ]s
esearch -db genome -query "${orgtaxsingle}"[TaxId] | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"

echo "Running esearch to save top filepath "

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




if [ -s "${orgnameidnospace}.genome.esearch.docsum" ];then


echo "${GENOME[@]}"

GENOME2=("${GENOME[0]}")

GENOME2=${GENOME2//ftp:/https:}

genomebase=(`basename $GENOME2`)
echo "${genomebase[0]}"
echo "${GENOME2[0]}"

wget -O "${orgnameidnospace}_genome.fna.gz" "${GENOME2[0]}"/"${genomebase[0]}"_genomic.fna.gz

echo " genome document summary ${GENOME2[0]}     ( ${orgnameid} ) "


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




fi
fi



if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 3 ]s

echo "No complete genomes identified in either the genomes or assembly databases, now trying any genome for                  ${orgnameid}"

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR genome[All Fields]) AND 2000:9999999[Slen]" | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"
echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 3 ]s

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR genome[All Fields]) AND 2000:9999999[Slen]" | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 1 ]s



echo " genome document summary  $GENOME                  ( ${orgnameid} ) "

echo " ${GENOME[@]} "

echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"




fi

if [ ! -s "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 2 )  + 3 ]s

echo "No RefSeq genome identified, now trying any genome                  ( ${orgnameid} ) "

sleep $[ ( $RANDOM % 2 )  + 3 ]s
esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR partial[All Fields]) AND 2000:9999999[Slen]" | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 3 ]s

esearch -db nuccore -query "${orgnameid}[All Fields] AND (complete[Text Word] OR partial[All Fields]) AND 2000:9999999[Slen]" | efetch -format docsum | tee "${orgnameidnospace}_wholedocsumdata.docsum"

echo " genome document summary $GENOME                  ( ${orgnameid} ) "

echo " ${GENOME[@]} "


echo ${orgtaxsingle} > "${orgnameidnospace}_TaxId.esearch.docsum"


fi



fi 




GENOME=()

echo " Completed search for viral organism number $c ${orgnameid}"
fi



if [ ! -s "${orgnameidnospace}_genome.fna.gz" ];then

rm "${orgnameidnospace}.genome.esearch.docsum"
rm "${orgnameidnospace}_genome.fna.gz"
rm "${orgnameidnospace}_TaxId.esearch.docsum"
rm "${orgnameidnospace}_wholedocsumdata.docsum"

fi

sleep $[ ( $RANDOM % 1 )  + 5 ]s

done

echo "All viral genomes available have now been downloaded"



touch ${programdir}${outdir}/finished.txt