#!/bin/bash
#SBATCH --account=OD-229285
#SBATCH --job-name filter_genomes
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=0:15:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1


while getopts 'i:o:p:z' c
do
  case $c in
    i) indir=$OPTARG ;;
    o) finishedsamp=$OPTARG ;;
    p) outdir=$OPTARG ;;
    z) samp=$OPTARG ;;
  esac
done

echo "output directory ${outdir} " 

echo "input directory ${indir} " 


# Have changed this to be a more broad and less strict filter system. Now I run it not against a genome but against sam hits 
# Threshold will be a wc -l with the setting being 100 reads

genome_hits=(`ls ${indir}/*hits_sam_no_header*`)

for seq in ${genome_hits[@]}; do
genomeref=(`awk '{sub(/.*\//, ""); sub(/_genome.*/, ""); sub(/.*_/, ""); print }'<<< $seq`)

echo "$genomeref"
echo " $samp "

echo "sequence $seq " && \
naligns=(`wc -l $seq`)
echo "number of reads aligned to reference genome $naligns "

if [ ${naligns} -ge 100 ] 
then
echo " detected significant hits for $genomeref. Now moving files to ${outdir} "
echo " taxid file to move ${indir}/*${genomeref}*TaxId* "
echo " bam file to move $seq "

cp $seq ${outdir}
cp ${indir}/*${genomeref}*TaxId* ${outdir}
fi
done


echo " finished copying"

touch ${finishedsamp}