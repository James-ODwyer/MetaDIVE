#!/bin/bash
#SBATCH --account=OD-229285 
#SBATCH --job-name snakemake_genome_assemble
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 12
#SBATCH --mem 8G
#SBATCH --time 24:00:00

eval "$(conda shell.bash hook)"

conda activate snakemake7

while getopts 'g:j:k:p:r:t:z:' c
do
  case $c in
    g) inputgenomedir=$OPTARG ;;
    j) RawSSUreads1=$OPTARG ;;
    k) RawSSUreads2=$OPTARG ;;
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

# I haven't updated the SSU reads name but they direct to earlier CO1 now (after I noticed some viral reads coming off at LSU). Doesn't change anything, just a note because I don't want to mess around with labelling 

# list of genomes in this sample
genomelist=(`ls ${programdir}${inputgenomedir}/*fna*`)

# get the correct sample raw reads to align
R1=(`ls ${programdir}${RawSSUreads1}`)
R2=(`ls ${programdir}${RawSSUreads2}`)

for fna in ${genomelist[@]}; do

genomeid=(`awk '{sub(/.*\//, ""); sub(/\.fna.gz/, ""); print }'<<< $fna`)

if [ ! -d ${programdir}${inputgenomedir}/${genomeid} ]; then mkdir -p ${programdir}${inputgenomedir}/${genomeid} 
fi
if [ ! -d ${programdir}${inputgenomedir}/${genomeid}/index ]; then mkdir -p ${programdir}${inputgenomedir}/${genomeid}/index 
fi
if [ ! -d ${programdir}${inputgenomedir}/${genomeid}/hits ]; then mkdir -p ${programdir}${inputgenomedir}/${genomeid}/hits 
fi
if [ ! -d ${programdir}${outputconsensusdir} ]; then mkdir -p ${programdir}${outputconsensusdir} 
fi

        bowtie2-build ${fna} ${programdir}${inputgenomedir}/${genomeid}/index/index_genome --threads ${SLURM_CPUS_PER_TASK} \
            --large-index && \
        bowtie2 -x ${programdir}${inputgenomedir}/${genomeid}/index/index_genome -1 ${R1} -2 ${R2} \
            --al-conc-gz ${programdir}${inputgenomedir}/${genomeid}/hits/${samplename}_${genomeid}_hits_R%.fastq.gz \
            -p 4 \
            -S ${programdir}${inputgenomedir}/${genomeid}/hits/${samplename}_${genomeid}_sam.sam \
            --very-sensitive-local
samtools view -F 4 -S ${programdir}${inputgenomedir}/${genomeid}/hits/${samplename}_${genomeid}_sam.sam > ${programdir}${inputgenomedir}/${samplename}_${genomeid}_hits_sam_no_header.bam

rm ${programdir}${inputgenomedir}/${genomeid}/hits/${samplename}_${genomeid}_sam.sam

done


touch ${programdir}${outputconsensusdir}/finished_align.txt
