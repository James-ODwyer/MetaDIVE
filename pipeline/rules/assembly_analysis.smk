# Only need to define wildcards agin if you want an if/else statement in defining them (May be useful for combined analysis. # will be useful for plans about host removal

def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]

assembly_to_run = config["Assembly_choice"]

Sensitivity = config["Sensitivity"]

def host_removed_dataR1(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R1_host_removed.fastq.gz"
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_unaligned.fastq.gz"
        ])

def host_removed_dataR2(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R2_host_removed.fastq.gz"
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_unaligned.fastq.gz"
        ])

def host_removed_dataR1_added_singles(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R1_host_removed_added_singles.fastq.gz",
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_unaligned_added_singles.fastq.gz"
        ])

def host_removed_dataR2_added_singles(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_R2_host_removed_added_singles.fastq.gz",
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_unaligned_added_singles.fastq.gz"
        ])

def host_removed_datasingles(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["host_remove_dir"] + "/{sample}_unaligned_singles.fastq.gz"
        ])
    if config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["SSU_dir"] + "/{sample}_unaligned_singles.fastq.gz"
        ])

def trinity_memtot(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "100"
        ])
    if config["Sensitivity"] == "High":
        return([
            "80"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "60"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "60"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "40"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "30"
        ])

def trinity_memtotMB(wildcards):
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


def trinity_mem_GCs(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "6"
        ])
    if config["Sensitivity"] == "High":
        return([
            "5"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "3"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "3"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "2"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "2"
        ])


def diamond_memtot(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "46"
        ])
    if config["Sensitivity"] == "High":
        return([
            "42"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "32"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "30"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "24"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "18"
        ])

def diamond_memtotMB(wildcards):
    if config["Sensitivity"] == "VHigh":
        return "46000"
    if config["Sensitivity"] == "High":
        return "42000"
    if config["Sensitivity"] == "Medium":
        return "32000"
    if config["Sensitivity"] == "Medium-Low":
        return "30000"
    if config["Sensitivity"] == "Low":
        return "24000"
    if config["Sensitivity"] == "Ultra_low":
        return "18000"


rule excluded_variable_host:
    input:
        R1 = host_removed_dataR1
    output:
        finishedfile = config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
    shell:
        """
        touch {output.finishedfile}
        """

rule excluded_variable:
    input:
        R1 = host_removed_dataR1
    output:
        finishedfile = "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
    shell:
        """
        touch {output.finishedfile}
        """

# Megahit is really poorly designed for snakemake. Snakemake creates all the directories around the program but megahit inputs a file directory
# And crashes if the directory exists still. You can get around this with --force but this command is hidden in the man page so don't forget etc

rule megahit_assembly:
    message:
        """
        Assembly step 1 
        Generating contigs from {wildcards.sample} using Megahit kmer based de novo assembly
        """
    input:
        R1 = host_removed_dataR1,
        R2 = host_removed_dataR2,
        singles = host_removed_datasingles
    output:
        contigsfile = config["sub_dirs"]["contig_dir_megahit"] + "/{sample}/final.contigs.fa"
    params:
        megahitpath = config["sub_dirs"]["contig_dir_megahit"] +"/{sample}",
        Assembly_size = config["Assembly_size"]
    log:
        "logs/" + config["sub_dirs"]["contig_dir_megahit"] + "/{sample}.log"
    threads: 16
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contig_dir_megahit"] + "/{sample}.txt"
    resources:
         mem_mb=32000
    shell:
        """
        megahit -1 {input.R1} -2 {input.R2} \
            --presets meta-sensitive \
            --memory 26000000 --force \
            -t {threads} \
            -o {params.megahitpath} \
            --min-contig-len {params.Assembly_size} \
            2> {log}
        """

# MEMDIR is of potential use here. If the trinity memory capacity was increased to the point of allowing for active storage of all files then this step could be run significantly faster
# Maybe as much as 2-5X given how many files the program actually generates
# Trinity has no option for incorporating paired and non paired reads without sacrificing read pair values (e.g., I can cat the reads all together and do as unaligned but give <0.5% of reads
# are typically singletons paired end information will be sacrified for 99.5% of reads for only small read number increases. ) I have opted to just not include the singletons into the trinity assembler
# Trinity is erroring out when casava format is retained. It appears to be a fastq format issue that trinity can't read the format properly. Easy fix is to add a /1 and /2 to the ends of every
# header line. 

rule trinity_assembly:
    message:
        """
        Assembly step 1 
        Generating contigs from {wildcards.sample} using Trinity kmer based de novo assembly
        """
    input:
        R1 = host_removed_dataR1,
        R2 = host_removed_dataR2
    output:
        contigsfiletemp = config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_trinity_temp.Trinity.fasta",
        R1_temp = temp(config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_R1_trinity.fastq"),
        R2_temp = temp(config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_R2_trinity.fastq"),
        contigsfile = config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_trinity.Trinity.fasta"
    params:
        temptrinitypath = config["Trinitytemppath"] + "/{sample}_trinity",
        trinity_file = ".Trinity.fasta",
        Assembly_size = config["Assembly_size"],
        trinitymemory = trinity_memtot,
        GCthreads = trinity_mem_GCs
    log:
        "logs/" + config["sub_dirs"]["contig_dir_trinity"] + "/{sample}.log"
    threads: 6
    priority: 10
    conda: "trinity"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contig_dir_trinity"] + "/{sample}.txt"
    resources:
        trinity=1,
        mem_mb=trinity_memtotMB(None)
    shell:
        """
        zcat {input.R1} | awk '/^@/ {{print $0"/1"}} !/^@/ {{print $0}}' > {output.R1_temp} && \
        zcat {input.R2} | awk '/^@/ {{print $0"/2"}} !/^@/ {{print $0}}' > {output.R2_temp} && \
        Trinity --left {output.R1_temp} --right {output.R2_temp} \
            --seqType fq \
            --full_cleanup \
            --CPU {threads} \
            --max_memory {params.trinitymemory}G \
            --min_contig_length {params.Assembly_size} \
            --normalize_reads \
            --bflyHeapSpaceMax 10G \
            --bflyGCThreads {params.GCthreads} \
            --output {params.temptrinitypath} \
            2> {log} && \
        mv {params.temptrinitypath}{params.trinity_file} {output.contigsfiletemp} && \
        reformat.sh in={output.contigsfiletemp} out={output.contigsfile} trd=t
        """

def Assembly_used(wildcards):
    if config["Assembly_choice"] == "Megahit":
        return([
            config["sub_dirs"]["contig_dir_megahit"] + "/{sample}/final.contigs.fa"
        ])
    if config["Assembly_choice"] == "Trinity":
        return([
            config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_trinity.Trinity.fasta"
        ])

rule Copy_contig_assemblies:
    input:
        Assemblysamples = Assembly_used
    output:
        contigsfile = config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
    shell:
        """
        cp {input.Assemblysamples} {output.contigsfile}
        """

def Assembly_used_initfilter(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_removed.fa"
        ])
    if config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
        ])

rule generate_kraken_contigs:
    message:
        """
        Generate additional viral assignments for contigs reads using kraken2
        for {wildcards.sample}
        """
    input:
        Assembly_used_initfilter
    output:
        krakenout = config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}_kraken_output.txt",
        krakenreport = config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}_kraken_report.txt",
        readslist2 = config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}_contigs_to_virus.txt"
    params:
        krakendb= config["Kraken_database"]
    log:
        "logs/" + config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}.log"
    threads: 2
    conda: "kraken2"
    resources:
        mem_mb=8000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}.txt"
    shell:
        """
        kraken2 --db {params.krakendb} \
            --threads {threads} \
            --minimum-hit-groups 7 \
            --report {output.krakenreport} \
            --output {output.krakenout} \
            {input} \
            2> {log} && \
        awk '$1 == "C" {{print $2}}' "{output.krakenout}" > {output.readslist2}
        """

# I had iterate on in the bash script but it appears the function was only introduced in 2.0.9. the latest anaconda version of Diamond is 2.0.8 
rule Diamond_blast:
    message:
        """
        Assembly step 2 
        Blast assembled contigs from {wildcards.sample} using Diamond blast
        """
    input:
        Assembly_used_initfilter
    output:
        matchesfile = config["sub_dirs"]["diamond_dir"] + "/{sample}_matches.m8",
        matchesfilefinished = config["sub_dirs"]["diamond_dir"] +"/{sample}_finished.txt"
    params:
        databasenr = config["diamond_database"],
        diamondmem = diamond_memtot
    log:
        "logs/" + config["sub_dirs"]["diamond_dir"] + "/{sample}.log"
    threads: 8
    priority: 10
    conda: "snakemake7"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["diamond_dir"] + "/{sample}.txt"
    resources:
        mem_mb=diamond_memtotMB(None)
    shell:
        """
        diamond blastx --db {params.databasenr} \
            --query {input} \
            --fast \
            --memory-limit {params.diamondmem} \
            -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp \
            --evalue 0.001 \
            --threads {threads} \
            -o {output.matchesfile} \
            --max-target-seqs 10 \
            2> {log} && \
        touch {output.matchesfile} && \
        touch {output.matchesfilefinished}
        """

rule analyse_diamond_hits:
    message:
        """
        Assembly step 3 
        generate taxonomies from each Diamond hit for {wildcards.sample} using R and taxonomizr package
        """
    input:
        diamondfile = config["sub_dirs"]["diamond_dir"] + "/{sample}_matches.m8",
        contigfile = Assembly_used,
        matchesfilefinished = config["sub_dirs"]["diamond_dir"] +"/{sample}_finished.txt"
    output:
        # change below to temp once completed everything
        contigsunassigned = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa",
        diamondfileupdated = config["sub_dirs"]["diamond_dir"] + "/{sample}_matches_updated.m8",
        contigsassignedtab = config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_Contigsallinformationassignment.txt",
        contigshostassigned = config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_host_aligned_contigs_list.txt",
        full_prot_host_id = config["sub_dirs"]["progress_main_pipeline_host"] + "/{sample}_host_id.txt"
    params:
        samplename = "{sample}",
        NCBI_nodes= config["NCBI_reference_dir"],
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["contigs_assigned"] + "/abundances/",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["contigs_assigned"],
        contigname = "{sample}_matches_unassigned.fa",
        full_prot_host_id_path = config["sub_dirs"]["progress_main_pipeline_host"]
    log:
        "logs/" + config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}.txt"
    threads: 1
    conda: "Rdataplotting"
    resources:
         mem_mb=3000
    priority: 10
    shell:
        """
        if [ ! -d {params.full_prot_host_id_path} ]; then mkdir -p {params.full_prot_host_id_path}; fi && \
        bash {config[program_dir]}scripts/get_taxIDs.sh -a {input.diamondfile} -i {params.NCBI_nodes} -o {output.diamondfileupdated} && \
        Rscript {config[program_dir]}scripts/extract_blastx_best_hits_HPC.R \
            --inputdiamond {output.diamondfileupdated} --inputcontig {input.contigfile} --name {params.samplename} --programdir {params.basedir} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.contigsunassigned} --abundances {params.abundances} --Log {log} \
            --savdir {params.wrkdir} --savcontig {params.contigname} --hosttaxid {output.full_prot_host_id}
        """

# Need to add a declare here to deal with the DNA vs not DNA assignment rules 
# added in the bam but Vamb recommends just filtering out errors and duplicates but keeping all unassigned reads. So that is what I have done (hence the weird sam flag 3584)

rule align_raws_to_contigs:
    message:
        """
        Assembly step 4 
        Aligning original reads to generated contigs for {wildcards.sample} using Bowtie2
        """
    input:
        Assembly = Assembly_used_initfilter,
        R1 = host_removed_dataR1,
        R2 = host_removed_dataR2,
        contigsassigned = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa"
    output:
        unalignedreads1 = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1.fastq.gz",
        unalignedreads2 = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2.fastq.gz",
        alignedreads1 = temp(config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1_aligned.fastq.gz"),
        alignedreads2 = temp(config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2_aligned.fastq.gz"),
        unalignedsingles = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutput = temp(config["sub_dirs"]["raws_to_contigs"] + "/{sample}.sam"),
        samoutputhits = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_hits.sam",
        samoutputhits2 = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_hits_3col.sam",
        bamoutput = config["sub_dirs"]["raws_to_contigs"] + "/{sample}.bam",
        samoutputunaligned = temp(config["sub_dirs"]["raws_to_contigs"] + "/{sample}_contigs_unaligned_reads.sam")
    params:
        indexpath = config["sub_dirs"]["temp_genome_path"] + "/{sample}/",
        temppathidx = config["sub_dirs"]["temp_genome_path"] + "/{sample}/{sample}_Tempidx",
        unalignedreads = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R%.fastq.gz",
        alignedreads = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["raws_to_contigs"] + "/{sample}.log"
    threads: 8
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_to_contigs"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.indexpath} ]; then mkdir -p {params.indexpath}; fi && \
        bowtie2-build {input.Assembly} {params.temppathidx} --threads {threads} \
            --large-index && \
        bowtie2 -x {params.temppathidx} -1 {input.R1} -2 {input.R2} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        samtools view -F 4 {output.samoutput} > {output.samoutputhits} && \
        cut -f 1-3 {output.samoutputhits} > {output.samoutputhits2} && \
        samtools view -b -F 3584 {output.samoutput} > {output.bamoutput} && \
        samtools view -@ {threads} -f 4 -h {output.samoutput} > {output.samoutputunaligned} && \
        samtools fastq -@ {threads} -1 {output.unalignedreads1} -2 {output.unalignedreads2} -s {output.unalignedsingles} -n {output.samoutputunaligned}
        """
