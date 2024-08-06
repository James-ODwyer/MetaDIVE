#!/bin/bash 
#SBATCH --account=OD-229285              # Required for CSIRO HPC. You need to specify your account. To see yours, write into a terminal (putty or other) get_project_codes
#SBATCH --job-name Build_databases_MetaDIVE      # named whatever you would like  
#SBATCH --nodes 1                        # nodes to use. 
#SBATCH --ntasks-per-node 1              # ntasks per node 
#SBATCH --cpus-per-task 4               # total number of CPUs to allocate.   
#SBATCH --mem 16G                       # Total memory. 
#SBATCH --time 6:00:00                 # Time requirements hh/mm/ss would recommend around 100 hours for large datasets. if it doesn't complete you can always launch the script again
#SBATCH --partition io


eval "$(conda shell.bash hook)"

conda activate snakemake7


databasedir=$(pwd)


mkdir -p bowtie2/LSU
mkdir -p bowtie2/SSU

mkdir LSU/


mv output_completed_taxonomyLSU.tsv.gz LSU/output_completed_taxonomyLSU.tsv.gz
mv LSU_combined.fasta.gz LSU/LSU_combined.fasta.gz

cd LSU/
cp LSU_combined.fasta.gz ../bowtie2/LSU/
gunzip output_completed_taxonomyLSU.tsv.gz
gunzip LSU_combined.fasta.gz


mmseqs createdb LSU_combined.fasta LSU_DB

# Download taxdump for creating taxonomy
mkdir ncbi-taxdump && cd ncbi-taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz
cd ..

# add taxonomy to LSU_DB
mmseqs createtaxdb LSU_DB tmp --ncbi-tax-dump ./ncbi-taxdump --tax-mapping-file output_completed_taxonomyLSU.tsv

mmseqs createindex LSU_DB tmp --search-type 3 --split 4

cd ..

mkdir SSU/


mv output_completed_taxonomySSU.tsv.gz SSU/output_completed_taxonomySSU.tsv.gz
mv SSU_combined.fasta.gz SSU/SSU_combined.fasta.gz



cd SSU/

cp SSU_combined.fasta.gz ../bowtie2/SSU/

gunzip output_completed_taxonomySSU.tsv.gz
gunzip SSU_combined.fasta.gz


mmseqs createdb SSU_combined.fasta SSU_DB

mkdir ncbi-taxdump && cd ncbi-taxdump

cp "$databasedir"/LSU/ncbi-taxdump/* .

cd ..

mmseqs createtaxdb SSU_DB tmp --ncbi-tax-dump ./ncbi-taxdump --tax-mapping-file output_completed_taxonomySSU.tsv


mmseqs createindex SSU_DB tmp --search-type 3 --split 6

cd ..
                

# Genomad
conda deactivate
conda activate genomad 

mkdir genomad

cd genomad

genomad download-database .
                   

cd ..
conda deactivate
conda activate snakemake7

# Bowtie indices 

cd bowtie2/LSU 


gunzip LSU_combined.fasta.gz


bowtie2-build "LSU_combined.fasta" LSU_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
            --large-index


cd ..
cd SSU


gunzip SSU_combined.fasta.gz


bowtie2-build "SSU_combined.fasta" SSU_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
            --large-index

cd ../..


# Create folder for accession taxonomy data (for Taxonkit plus taxonomizr)

mkdir taxonomy

mv nameNode.sqlite.gz taxonomy/

cd taxonomy

gunzip nameNode.sqlite.gz

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O taxdump.tar.gz

tar -xzvf taxdump.tar.gz names.dmp nodes.dmp


cd ..


# Create folder for adapter sequences (provided as gzip file)
mkdir adapters

mv trueseq_PE_adaptors_shortlist.fa.gz adapters/

cd adapters

gunzip trueseq_PE_adaptors_shortlist.fa.gz

cd ..





mkdir Phix_reference

mv PhiX_genome.fa.gz Phix_reference/

cd Phix_reference 

gunzip PhiX_genome.fa.gz

bowtie2-build "PhiX_genome.fa" Phix_reference_idx --threads ${SLURM_CPUS_PER_TASK}

cd ..



# CO1 reference genome generation for MMseqs2

conda activate snakemake7

mkdir CO1
mkdir -p bowtie2/CO1

mv output_completed_taxonomyCO1.tsv.gz CO1/output_completed_taxonomyCO1.tsv.gz
mv CO1_all_references.fasta.gz CO1/CO1_all_references.fasta.gz

cd CO1/
cp CO1_all_references.fasta.gz ../bowtie2/CO1/
gunzip output_completed_taxonomyCO1.tsv.gz
gunzip CO1_all_references.fasta.gz


mmseqs createdb CO1_all_references.fasta CO1_DB

# Download taxdump for creating taxonomy
mkdir ncbi-taxdump && cd ncbi-taxdump

cp "$databasedir"/LSU/ncbi-taxdump/* .

cd ..

# add taxonomy to CO1_DB
mmseqs createtaxdb CO1_DB tmp --ncbi-tax-dump ./ncbi-taxdump --tax-mapping-file output_completed_taxonomyCO1.tsv


mmseqs createindex CO1_DB tmp --search-type 3 --split 4

cd ..



# CO1 bowtie2


cd bowtie2/CO1

gunzip CO1_all_references.fasta.gz

bowtie2-build "CO1_all_references.fasta" CO1_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
            --large-index

cd ..
cd ..


conda activate Rdataplotting

Rscript install_d3Tree.R

conda deactivate





