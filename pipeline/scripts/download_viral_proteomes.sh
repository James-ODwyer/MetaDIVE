#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_genome_download
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 4:00:00
#SBATCH --partition io

eval "$(conda shell.bash hook)"

conda activate snakemake7


# read in required inputs
# input is the annotated gene file from Dramv and Matt's py code to combine the annot file with the gff gene file
# outdir will be the 16_REFERENCE_VIRAL_GENOMES/{sample}/
# programdir will be the config["programdir"]


############### NOTE!!!!!!!!!!!!!! The tsv files from R can be new line from Windows format. The files were all failing because of this newline stopping linux from correctly calling the 
# variables # just preemptively fix from now on will probably be easier 


while getopts 'i:p:z:' c
do
  case $c in
    i) input=$OPTARG ;;
    p) outdir=$OPTARG ;;
    z) programdir=$OPTARG ;;
  esac
done


dos2unix $input

# api key here, but should make it a variable in the config file and feed it in in future.
export NCBI_API_KEY=b81d4cea18edf2ff96dd6aeb096d6eb56009


declare -a orgname=()
declare -a orgnamenospace=()

IFS=$'\n'
while read line; do 
echo " full line ${line[@]} "
orgname1=(`echo "${line[@]}"`)
echo " orgname 1 ${orgname1} "
orgname1a="${orgname1[@]}"
orgname2=(${orgname1a[@]// /})
echo " ${orgname2} "
orgname+=(" ${orgname1a} ")
orgnamenospace+=("${orgname2}")
done < sampleX_viralspecies_and_secondary.tsv



START=0
END=($((${#orgnamenospace[@]}-1)))

echo $START
echo $END

for (( c=$START; c<=$END; c++ ))
do

orgnameid=(`echo "${orgname[$c]}"`)
echo " $orgnameid "
orgnameidnospace=(`echo -n "${orgnameid//[[:space:]]/}"`)

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 1 )  + 1 ]s

# Step 1 is there a full genome

echo " Running esearch to generate docsum  for ${orgnameid} "
esearch -db assembly -query "${orgnameid}"[ORGN] -spell | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | tee "${orgnameidnospace}.genome.esearch.docsum"

esearch -db assembly -query "${orgnameid}"[ORGN] -spell | efetch -rettype fasta_cds_translate -retmode text
echo " Running esearch to save top filepath "

#esearch -db assembly -query "${orgnamesp1[$c]}"[ORGN] | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | tee "${orgnameidnospace}.genome.esearch.docsum"

[ -s "${orgnameidnospace}.genome.esearch.docsum" ] || rm "${orgnameidnospace}.genome.esearch.docsum"

fi


if [ -f "${orgnameidnospace}.genome.esearch.docsum" ];then

echo " no refseq genome assembly avaliable, will search nuccore for all nucleotide segments for ${orgnameid} "

sleep $[ ( $RANDOM % 1 )  + 3 ]s

esearch -db assembly -query "${orgnameid}"[ORGN] | efetch -format docsum | tee "${orgnameidnospace}.genome.esearch.docsum"

[ -s "${orgnameidnospace}.genome.esearch.docsum" ] || rm "${orgnameidnospace}.genome.esearch.docsum"

fi


if [ -f "${orgnameidnospace}.genome.esearch.docsum" ];then

echo " No genome was found for ${orgnameid} will now download information on avaliable nucleotide sequences "

sleep $[ ( $RANDOM % 1 )  + 2 ]s

esearch -db nuccore -query "${orgnameid}"[ORGN] | efetch -format docsum | tee "${orgnameidnospace}.genome.esearch.docsum"

[ -s ${orgnameidnospace}.genome.esearch.docsum ] || rm ${orgnameidnospace}.genome.esearch.docsum

fi

done
# something like this! 
esearch -db nuccore -query "${orgnamesp1[$c]}" | efetch  efetch -format gb > "${orgnameidnospace}_genome.fna.gz"

#fasta_cds_aa

sleep $[ ( $RANDOM % 2 )  + 2 ]s

GENOME=$(esearch -db assembly -query "${orgnameid}"[ORGN] | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | \
awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}')


echo " genome document summary  $GENOME "

echo " ${GENOME[@]} "

fi

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 5 )  + 2 ]s

echo "No RefSeq genome identified, now trying any genome"
esearch -db genome -query "${orgnameid}"[ORGN] | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | tee "${orgnameidnospace}.genome.esearch.docsum"
echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 2 ]s

GENOME=$(esearch -db assembly -query "${orgnameid}"[ORGN] | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank  | \
awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}')


echo " genome document summary  $GENOME "

echo " ${GENOME[@]} "

fi

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 5 )  + 2 ]s

echo "No RefSeq genome identified, now trying any genome"
GENOME=$(esearch -db nuccore -query '("${orgnameid}"[Organism] OR "${orgnameid}"[All Fields]) AND complete cds[Text Word]' | efetch -format docsum | tee "${orgnameidnospace}.genome.esearch.docsum")
echo "Running esearch to save top filepath "

sleep $[ ( $RANDOM % 2 )  + 2 ]s

esearch -db nuccore -query '("${orgnameid}"[Organism] OR "${orgnameid}"[All Fields]) AND complete cds[Text Word]' | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 2 ]s



echo " genome document summary  $GENOME "

echo " ${GENOME[@]} "

fi

if [ ! -f "${orgnameidnospace}.genome.esearch.docsum" ];then

sleep $[ ( $RANDOM % 5 )  + 2 ]s

echo "No RefSeq genome identified, now trying any genome"
GENOME=$(esearch -db nuccore -query '("${orgnameid}"[Organism] OR "${orgnameid}"[All Fields]) AND partial genome[Text Word]' | efetch -format docsum | tee "${orgnameidnospace}.genome.esearch.docsum")

sleep $[ ( $RANDOM % 2 )  + 2 ]s
esearch -db nuccore -query '("${orgnameid}"[Organism] OR "${orgnameid}"[All Fields]) AND "partial genome"[Text Word]' | efetch -format fasta > "${orgnameidnospace}_genome.fna.gz"

sleep $[ ( $RANDOM % 2 )  + 2 ]s


echo " genome document summary $GENOME "

echo " ${GENOME[@]} "

fi



if [ ${#GENOME} -ge 2 ]; then

if [ ! -f "${orgnameidnospace}_genome.fna.gz" ]; then

sleep $[ ( $RANDOM % 10 )  + 2 ]s

echo " downloading ${GENOME[0]} "

wget ${GENOME[0]}

file=(`awk '{sub(/.*\//, ""); print }'<<< ${GENOME[0]}`)

mv ${file} "${orgnameidnospace}_genome.fna.gz"

fi

fi

GENOME=()

echo " Completed search for viral gene number $c ${orgnameid[$c]}"


done

echo "All viral genomes available have now been downloaded"