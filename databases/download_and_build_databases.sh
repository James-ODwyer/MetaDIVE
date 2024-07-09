#!/bin/bash 
#SBATCH --account=OD-229285              # Required for CSIRO HPC. You need to specify your account. To see yours, write into a terminal (putty or other) get_project_codes
#SBATCH --job-name Build_databases_MetaDIVE      # named whatever you would like  
#SBATCH --nodes 1                        # nodes to use. 
#SBATCH --ntasks-per-node 1              # ntasks per node 
#SBATCH --cpus-per-task 4               # total number of CPUs to allocate.   
#SBATCH --mem 16G                       # Total memory. 
#SBATCH --time 8:00:00                 # Time requirements hh/mm/ss would recommend around 100 hours for large datasets. if it doesn't complete you can always launch the script again
#SBATCH --partition io


eval "$(conda shell.bash hook)"

conda activate snakemake7

mkdir -p bowtie2/LSU
mkdir -p bowtie2/SSU

mkdir LSU/

cd LSU/

# build name.dmp, node.dmp from SILVA taxonomy
mkdir taxonomy/ && cd "$_"
wget -O SILVA_LSU_NR99_taxmap_138_1.txt.gz "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz"
buildNCBITax=$(cat << 'EOF'
BEGIN{
ids["root"]=1;
print "1\t|\t1\t|\tno rank\t|\t-\t|" > "nodes.dmp"
print "1\t|\troot\t|\t-\t|\tscientific name\t|" > "names.dmp";
}
{ n=split($1, a, ";");
gsub("domain", "superkingdom", $3);
ids[$1]=$2;
gsub(/[^,;]*;$/,"",$1);
pname=$1;
if(n==2){
pname="root"
}
pid=ids[pname];
printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", $2, pid, $3) > "nodes.dmp";
printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n",$2,a[n-1]) >"names.dmp";
}
EOF
)
awk -F'\t' "$buildNCBITax" <(gunzip -c SILVA_LSU_NR99_taxmap_138_1.txt.gz)
touch merged.dmp
touch delnodes.dmp
cd ..
# create the database SILVA database from Nr99 fasta
wget -O SILVA_LSU_NR99_138_1.gz "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
cp SILVA_LSU_NR99_138_1.gz ../bowtie2/LSU/
mmseqs createdb SILVA_LSU_NR99_138_1.gz SILVA_DB
# add taxonomy to SILVA_DB
wget -O SILVA_LSU_NR99_taxmap_138_1.acc_taxid.gz "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_lsu_138.1.acc_taxid.gz"
gunzip SILVA_LSU_NR99_taxmap_138_1.acc_taxid.gz
mmseqs createtaxdb SILVA_DB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file SILVA_LSU_NR99_taxmap_138_1.acc_taxid

mmseqs createindex SILVA_DB tmp --search-type 3 --split 4

cd ..

mkdir SSU/

cd SSU/

# build name.dmp, node.dmp from SILVA taxonomy
mkdir taxonomy/ && cd "$_"
wget -O SILVA_SSU_NR99_taxmap_138_1.txt.gz "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz"
buildNCBITax=$(cat << 'EOF'
BEGIN{
ids["root"]=1;
print "1\t|\t1\t|\tno rank\t|\t-\t|" > "nodes.dmp"
print "1\t|\troot\t|\t-\t|\tscientific name\t|" > "names.dmp";
}
{ n=split($1, a, ";");
gsub("domain", "superkingdom", $3);
ids[$1]=$2;
gsub(/[^,;]*;$/,"",$1);
pname=$1;
if(n==2){
pname="root"
}
pid=ids[pname];
printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", $2, pid, $3) > "nodes.dmp";
printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n",$2,a[n-1]) >"names.dmp";
}
EOF
)
awk -F'\t' "$buildNCBITax" <(gunzip -c SILVA_SSU_NR99_taxmap_138_1.txt.gz)
touch merged.dmp
touch delnodes.dmp
cd ..
# create the database SILVA database from Nr99 fasta
wget -O SILVA_SSU_NR99_138_1.gz "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
cp SILVA_SSU_NR99_138_1.gz ../bowtie2/SSU/
mmseqs createdb SILVA_SSU_NR99_138_1.gz SILVA_DB
# add taxonomy to SILVA_DB
wget -O SILVA_SSU_NR99_taxmap_138_1.acc_taxid.gz "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.1.acc_taxid.gz"
gunzip SILVA_SSU_NR99_taxmap_138_1.acc_taxid.gz
mmseqs createtaxdb SILVA_DB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file SILVA_SSU_NR99_taxmap_138_1.acc_taxid

mmseqs createindex SILVA_DB tmp --search-type 3 --split 4

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

gunzip SILVA_LSU_NR99_138_1.gz


bowtie2-build "SILVA_LSU_NR99_138_1" LSU_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
            --large-index


cd ..
cd SSU

gunzip SILVA_SSU_NR99_138_1.gz


bowtie2-build "SILVA_SSU_NR99_138_1" SSU_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
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
mkdir CO1/taxonomy
cp CO1_reference_sequences_Metazoa.fasta.gz bowtie2/CO1/
cp tax-mapping-file accession_taxid_mapping.tsv.gz

gunzip CO1_reference_sequences_Metazoa.fasta.gz
gunzip tax-mapping-file accession_taxid_mapping.tsv.gz
mv CO1_reference_sequences_Metazoa.fasta CO1/

cp taxonomy/names.dmp CO1/taxonomy/names.dmp

cp taxonomy/nodes.dmp CO1/taxonomy/nodes.dmp

touch CO1/taxonomy/merged.dmp
touch CO1/taxonomy/delnodes.dmp

cd CO1

mmseqs createdb CO1_reference_sequences_Metazoa.fasta CO1_DB

mmseqs createtaxdb CO1_DB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file accession_taxid_mapping.tsv

mmseqs createindex CO1_DB tmp --search-type 3 --split 4

cd ..




# CO1 bowtie2



cd bowtie2/CO1

gunzip CO1_reference_sequences_Metazoa.fasta.gz


bowtie2-build "CO1_reference_sequences_Metazoa.fasta" CO1_reference_idx --threads ${SLURM_CPUS_PER_TASK} \
            --large-index

cd ..
cd ..






