def host_removed_dataR1(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R1_host_removed.fastq.gz"
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R1.fastq.gz"
        ])

def host_removed_dataR2(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R2_host_removed.fastq.gz"
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R2.fastq.gz"
        ])

def Assembly_used_initfilter(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_removed.fa"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
        ])

def hostdetected(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/{sample}_top_host_species_overall.txt"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])

def blastnreturn(wildcards):
    if config["DNA_assign_blastn"] == 'yes':
        return([
            config["sub_dirs"]["contigs_assigned_nucl"] + "/abundances/{sample}_host_aligned_contigs_list.txt"
        ])
    elif config["DNA_assign_blastn"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])


def genomad_mem_tot(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "40"
        ])
    if config["Sensitivity"] == "High":
        return([
            "36"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "32"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "26"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "24"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "20"
        ])

def genomad_memtotMB(wildcards):
    if config["Sensitivity"] == "VHigh":
        return "100000"
    if config["Sensitivity"] == "High":
        return "80000"
    if config["Sensitivity"] == "Medium":
        return "60000"
    if config["Sensitivity"] == "Medium-Low":
        return "60000"
    if config["Sensitivity"] == "Low":
        return "40000"
    if config["Sensitivity"] == "Ultra_low":
        return "30000"



def genomad_mem_splits(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "18"
        ])
    if config["Sensitivity"] == "High":
        return([
            "28"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "32"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "42"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "64"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "76"
        ])



rule subset_contig_files:
    message:
        """
        Genomad step 0A  
        Subsetting contigs file to remove contigs assigned as the host genus already for {wildcards.sample}
        """
    input:
        assembly_results = Assembly_used_initfilter,
        blastnrun = blastnreturn,
        contigshostassigned = config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_host_aligned_contigs_list.txt"
    output:
        contigs_nohost = config["sub_dirs"]["Genomad_viral_ids_prep"] + "/{sample}_contigs.fasta"
    params:
        contiglist = config["sub_dirs"]["Genomad_viral_ids_prep"] + "/{sample}_combined_host_contigs_list.txt"
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_ids_prep"] + "/{sample}.log"
    threads: 2
    priority: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}.txt"
    resources:
        mem_mb=2000
    shell:
        """
        if ! [ -s {input.blastnrun} ];then
        cat {input.contigshostassigned} > {params.contiglist}
        fi && \
        if [ -s {input.blastnrun} ];then
        cat {input.contigshostassigned} {input.blastnrun} > {params.contiglist}
        fi && \
        seqkit grep -vf {params.contiglist} {input.assembly_results} > {output.contigs_nohost}
        """

rule genomad_contigs:
    message:
        """
        Genomad step 1 
        generating probable viral contigs using neural network similarities for {wildcards.sample}
        """
    input:
        contigs_nohost = config["sub_dirs"]["Genomad_viral_ids_prep"] + "/{sample}_contigs.fasta"
    output:
        outfiletemp=temp(config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}_contigs.fasta"),
        viral_contigsp = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_proteins.faa",
        viral_contigsn = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna",
        FDRrates = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv"
    params:
        outdir = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}",
        databasegenomad = config["Genomaddb"],
        genomad_mem_splitnum = genomad_mem_splits,
        genomad_mem = genomad_mem_tot
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}.log"
    threads: 8
    priority: 1
    conda: "genomad"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}.txt"
    resources:
        mem_mb=genomad_memtotMB(None)
    shell:
        """
        cp {input.contigs_nohost} {output.outfiletemp} && \
        genomad end-to-end --min-score 0.75 --max-fdr 0.1 --min-plasmid-hallmarks-short-seqs 0 --min-virus-hallmarks-short-seqs 0 --max-uscg 50 --enable-score-calibration --splits {params.genomad_mem_splitnum} --cleanup {output.outfiletemp} {params.outdir} {params.databasegenomad}
        """

rule diamondblastn_viralhits:
    message:
        """
        Genomad step 2 
        Diamond blastn identified genomad viral contigs from {wildcards.sample} using blastn
        """
    input:
        viral_contigs = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna"
    output:
        blastfile = config["sub_dirs"]["Genomad_viral_blastn"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["Genomad_viral_blastn"]
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_blastn"] + "/{sample}.log"
    threads: 6
    priority: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_blastn"] + "/{sample}.txt"
    resources:
        mem_mb=20000
    shell:
        """
        blastn -query {input.viral_contigs} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 1 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -word_size 22 \
            -out {output.blastfile} \
            2> {log}
        touch {output.blastfile}
        """

rule remove_host_hits_blastn:
    message:
        """
        Genomad step 3 
        generate taxonomies from each Diamond hit for {wildcards.sample} using R and taxonomizr package and identify host associated contigs
        """
    input:
        hostsp = hostdetected,
        blastfile = config["sub_dirs"]["Genomad_viral_blastn"] + "/{sample}_matches_nucleotide.m8",
        full_prot_host_id = config["sub_dirs"]["progress_main_pipeline_host"] + "/{sample}_host_id.txt",
        viral_contigs = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna"
    output:
        hostcontigs = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_host_aligned_contigs.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["Genomad_viral_host_rem"]
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}.log"
    threads: 1
    conda: "Rdataplotting"
    priority: 1
    shell:
        """
        num_rows=$(wc -l < "{input.blastfile}") && \
        if [ "$num_rows" -gt 1 ]; then
        Rscript {config[program_dir]}scripts/Identify_host_contigs_genomad_HPC.R \
            --inputdiamond {input.blastfile} --name {params.samplename} --programdir {params.basedir} \
            --Accnode {params.namenodedatabase} --output {output.hostcontigs} --Log {log} \
            --savdir {params.wrkdir} --hostsp {input.hostsp} --contighostsp {input.full_prot_host_id}
        fi && \
        touch {output.hostcontigs}
        """

#Note if this doesn't work because there are non hits in the bam file, add an additional samtools view to get to just hits.

rule generate_depth_bams_for_genome_binning:
    message:
        """
        Genomad step 6 
        aligning raw reads to viral contigs from genomad and determine depth statistics for {wildcards.sample} using bowtie2 and metabat2
        """
    input:
        viral_contigs = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna",
        R1 = host_removed_dataR1,
        R2 = host_removed_dataR2,
        contigsassigned = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa",
        hostcontigs = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_host_aligned_contigs.txt"
    output:
        contigsnohost = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_contigs_virus_only.fna",
        alignedreads1 = temp(config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_R1_aligned.fastq.gz"),
        alignedreads2 = temp(config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_R2_aligned.fastq.gz"),
        samoutput = temp(config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}.sam"),
        bamoutputsorted = config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_hits_sorted.bam",
        metabat_depths = config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_metabat_depth_stats.txt"
    params:
        indexpath = config["sub_dirs"]["temp_genome_path"] + "/{sample}_VIRAL_CONTIGS/",
        temppathidx = config["sub_dirs"]["temp_genome_path"] + "/{sample}_VIRAL_CONTIGS/{sample}_Tempidx",
        unalignedreads = config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_R%.fastq.gz",
        alignedreads = config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}.log"
    threads: 6
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.indexpath} ]; then mkdir -p {params.indexpath}; fi && \
        seqkit grep -vf {input.hostcontigs} {input.viral_contigs} > {output.contigsnohost} && \
        bowtie2-build {output.contigsnohost} {params.temppathidx} --threads {threads} \
            --large-index && \
        bowtie2 -x {params.temppathidx} -1 {input.R1} -2 {input.R2} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --sensitive \
            2> {log} && \
        samtools sort -@ {threads} -o {output.bamoutputsorted} {output.samoutput} && \
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat_depths} {output.bamoutputsorted}
        """

rule Metabat_binning:
    message:
        """
        Genomad step 7 
        genome binning of viral like contigs for {wildcards.sample} and Metabat2
        """
    input:
        contigsnohost = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_contigs_virus_only.fna",
        metabat_depths = config["sub_dirs"]["Genomad_viral_binning_prep"] + "/{sample}_metabat_depth_stats.txt"
    output:
        outfile = config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}/finished.txt"
    params:
        outpath = config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}/{sample}"
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}.log"
    threads: 8
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}.txt"
    shell:
        """
        metabat2 -t {threads} -m 1500 --minS 50 --maxEdges 400 --minClsSize 1600 \
            -i {input.contigsnohost} \
            -a {input.metabat_depths} \
            --unbinned \
            -o {params.outpath} \
           2> {log} && \
        touch {output.outfile}
        """

rule sort_metabat_results:
    message:
        """
        Genomad step 8 
        Sort the completed Metabat files, and naming all contigs through genomad or diamond classifications
        """
    input:
        metabat_finished = config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}/finished.txt",
        summaryrawdiamondhitsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt",
        FDRrates = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv"
    output:
        outfile = config["sub_dirs"]["Genomad_metabat_sorted"] + "/{sample}/finished.txt"
    params:
        inpath = config["sub_dirs"]["Genomad_viral_binning_metabat"] + "/{sample}",
        outpath = config["sub_dirs"]["Genomad_metabat_sorted"] + "/{sample}"
    log:
        "logs/" + config["sub_dirs"]["Genomad_metabat_sorted"] + "/{sample}.log"
    threads: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_metabat_sorted"] + "/{sample}.txt"
    shell:
        """
        bash {config[program_dir]}scripts/sort_metabat_sequences.sh -j {input.FDRrates} -p {input.summaryrawdiamondhitsfiles} -s {params.inpath} -z {params.outpath} && \
        touch {output.outfile}
        """

rule blastp_viral_hits:
    message:
        """
        Genomad step 4 
        Diamond blastp identified genomad viral contigs from {wildcards.sample} using Diamond blast
        """
    input:
        viral_contigs = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_proteins.faa",
        nonhostcontigs = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_non_host_aligned_contigs.tsv"
    output:
        subset_contigs = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_contigs_virus_proteins_non_host.faa",
        matchesfile = config["sub_dirs"]["Genomad_viral_blastp"] + "/{sample}_matches.m8"
    params:
        databasevir = config["virus_protein_database"]
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_blastp"] + "/{sample}.log"
    threads: 16
    priority: 1
    conda: "snakemake7"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Genomad_viral_blastp"] + "/{sample}.txt"
    resources:
        mem_mb=16000
    shell:
        """
        seqtk subseq {input.viral_contigs} {input.nonhostcontigs} > {output.subset_contigs} && \
        blastp -db {params.databasevir} \
            -query {output.subset_contigs} \
            -evalue 0.001 -culling_limit 3 -max_hsps 1 \
            -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend sstart send staxids stitle qcovhsp' \
            -num_threads {threads} \
            -out {output.matchesfile}
        """

rule remove_host_hits:
    message:
        """
        Genomad step 3 
        generate taxonomies from each Diamond hit for {wildcards.sample} using R and taxonomizr package and identify host associated contigs
        """
    input:
        hostsp = hostdetected,
        diamondfile = config["sub_dirs"]["Genomad_viral_blastp"] + "/{sample}_matches.m8",
    output:
        nonhostcontigs = config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}_non_host_aligned_contigs.tsv"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["Genomad_viral_host_rem"]
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_host_rem"] + "/{sample}.log"
    threads: 1
    conda: "Rdataplotting"
    priority: 1
    shell:
        """
        Rscript {config[program_dir]}scripts/Identify_host_contigs_genomad_HPC.R \
            --inputdiamond {input.diamondfile} --name {params.samplename} --programdir {params.basedir} \
            --Accnode {params.namenodedatabase} --output {output.nonhostcontigs} --Log {log} \
            --savdir {params.wrkdir} --hostsp {input.hostsp}
        """


rule generate_top_viral_lists:
    message:
        """
        Genomad step 5 
        generate combined gene and contig lists of different qualities from Genomad for {wildcards.sample}
        """
    input:
        diamondfile = config["sub_dirs"]["Genomad_viral_blastp"] + "/{sample}_matches.m8",
        FDRrates = config["sub_dirs"]["Genomad_viral_ids"] + "/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv"
    output:
        summary_viral_hits = config["sub_dirs"]["Genomad_viral_categorise_species"] + "/{sample}_summary_all_viral_species.tsv"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        namenodedatabase = config["Accession_allnamenode"],
        wrkdir = config["sub_dirs"]["Genomad_viral_categorise_species"]
    log:
        "logs/" + config["sub_dirs"]["Genomad_viral_categorise_species"] + "/{sample}.log"
    threads: 2
    conda: "Rdataplotting"
    priority: 1
    shell:
        """
        length=(`wc -l {input.diamondfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_viral_species_multi_genes_genomad.R \
            --inputdiamond {input.diamondfile} --name {params.samplename} --programdir {params.basedir} \
            --output {output.summary_viral_hits} --Log {log} \
            --savdir {params.wrkdir} --FDRrates {input.FDRrates} --Accnode {params.namenodedatabase}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        touch {output.summary_viral_hits}
        fi
        """

