
def raws_after_hostR1(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R1.fastq.gz"
        ])
    elif config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1.fastq.gz"
        ])


def raws_after_hostR2(wildcards):
    if config["Host_filter"] == "yes":
        return([
            config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R2.fastq.gz"
        ])
    elif config["Host_filter"] == "no":
        return([
            config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2.fastq.gz"
        ])

def Diamond_memraw(wildcards):
    if config["Sensitivity"] == "VHigh":
        return([
            "36"
        ])
    if config["Sensitivity"] == "High":
        return([
            "34"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "28"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "24"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "20"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "16"
        ])


def Diamond_memrawMB(wildcards):
    if config["Sensitivity"] == "VHigh":
        return "36000"
    if config["Sensitivity"] == "High":
        return "34000"
    if config["Sensitivity"] == "Medium":
        return "28000"
    if config["Sensitivity"] == "Medium-Low":
        return "24000"
    if config["Sensitivity"] == "Low":
        return "20000"
    if config["Sensitivity"] == "Ultra_low":
        return "16000"

rule unmapped_reads_diamond:
    message:
        """
        Assembly step 5 
        Blast unaligned raw reads from {wildcards.sample} using Diamond blast
        """
    input:
        unalignedreads1 = raws_after_hostR1,
        unalignedreads2 = raws_after_hostR2
    output:
        interleavedfiles = temp(config["sub_dirs"]["raws_to_contigs"] + "/{sample}_interleaved.fastq.gz"),
        diamondfile = config["sub_dirs"]["diamond_raws"] + "/{sample}_matches.m8"
    params:
        databasenr = config["diamond_database"],
        unalignedreads1_prior = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1.fastq.gz",
        unalignedreads2_prior = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2.fastq.gz",
        diamondmem = Diamond_memraw
    log:
        "logs/" + config["sub_dirs"]["diamond_raws"] + "/{sample}.txt"
    threads: 8
    priority: 5
    conda: "snakemake7"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["diamond_raws"] + "/{sample}.txt"
    resources:
        mem_mb=Diamond_memrawMB(None)
    shell:
        """
        length=(`wc -l {input.unalignedreads1}`) && \
        if [ ${{length}} -ge 1 ]
        then        
        seqfu ilv -1 {input.unalignedreads1} -2 {input.unalignedreads2} -o {output.interleavedfiles} && \
        diamond blastx --db {params.databasenr} \
            --query {output.interleavedfiles} \
            --fast \
            --max-target-seqs 1 \
            -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp \
            --evalue 0.0001 \
            --threads {threads} \
            -o {output.diamondfile} \
            --memory-limit {params.diamondmem} \
            --masking 0 \
            2> {log}
        fi  && \
        if [ ${{length}} -eq 0 ]
        then
        seqfu ilv -1 {params.unalignedreads1_prior} -2 {params.unalignedreads2_prior} -o {output.interleavedfiles} && \
        diamond blastx --db {params.databasenr} \
            --query {output.interleavedfiles} \
            --fast \
            --max-target-seqs 1 \
            -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp \
            --evalue 0.0001 \
            --threads {threads} \
            -o {output.diamondfile} \
            --memory-limit {params.diamondmem} \
            --masking 0 \
            2> {log}
        fi        
        """

rule analyse_diamond_hits_raws:
    message:
        """
        Assembly step 6 
        generate taxonomies from each Diamond hit for {wildcards.sample} raw reads using R and taxonomizr package
        """
    input:
        diamondfile = config["sub_dirs"]["diamond_raws"] + "/{sample}_matches.m8"
    output:
        diamondfileupdated = config["sub_dirs"]["diamond_raws"] + "/{sample}_matches_updated.m8",
        sampfinished = config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}_finished.txt",
        superkingdoms = config["sub_dirs"]["contigs_assigned_raw"] + "/abundances/{sample}_superkingdoms.txt",
        contigsassigned = config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}_contigs_file.txt"
    params:
        samplename = "{sample}",
        NCBI_nodes= config["NCBI_reference_dir"],
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["contigs_assigned_raw"] + "/abundances/",
        contigsassigned = config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}_contigs_file.txt"
    log:
        "logs/" + config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    priority: 20
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}.txt"
    resources:
        mem_mb=6000
    shell:
        """
        bash {config[program_dir]}scripts/get_taxIDs.sh -a {input.diamondfile} -i {params.NCBI_nodes} -o {output.diamondfileupdated} && \
        Rscript {config[program_dir]}scripts/extract_blastx_best_hits_HPC_raws.R \
            --inputdiamond {output.diamondfileupdated} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.contigsassigned} --abundances {params.abundances} && \
        touch {output.sampfinished} && \
        touch {output.superkingdoms} && \
        touch {output.contigsassigned}
        """

rule check_results_in_blastn:
    message:
        """
        Running raw reads assigned to viruses in blast x against blastn database as false positive check.
        Blasting reads of {wildcards.sample} (blastn)
        """
    input:
        R1 = raws_after_hostR1,
        R2 = raws_after_hostR2,
        readslist = config["sub_dirs"]["Summary_results"] + "/{sample}_raw_read_names_to_virus.txt"
    output:
        fastqreads = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads.fastq",
        fastareads = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads.fasta",
        fasta_clust = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads_clust.fasta",
        fasta_clust_ref_file = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads_clust.fasta.clstr",
        blastfile = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["contigs_assigned_nucl"]
    log:
        "logs/" + config["sub_dirs"]["raws_blastn_check"] + "/{sample}_blastn.log"
    threads: 4
    resources:
        mem_mb=16000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_blastn_check"] + "/{sample}_blastn.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        lengthunassigned=(`wc -l {input.readslist}`) && \
        if [ ${{lengthunassigned}} -ge 1 ]
        then
        zcat {input.R1} | grep --no-group-separator -A 3 -F -f "{input.readslist}" > {output.fastqreads}
        zcat {input.R2} | grep --no-group-separator -A 3 -F -f "{input.readslist}" >> {output.fastqreads}
        seqtk seq -a {output.fastqreads} > {output.fastareads}
        cd-hit -i {output.fastareads} -o {output.fasta_clust} -c 0.95 -n 5 -T {threads} -d 0 -M 14000
        blastn -query {output.fasta_clust} \
            -db {params.blastdb} \
            -evalue 0.01 \
            -max_target_seqs 1 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -word_size 21 \
            -out {output.blastfile} \
            2> {log}
        fi && \
        touch {output.blastfile} && \
        touch {output.fastqreads} && \
        touch {output.fastareads} && \
        touch {output.fasta_clust} && \
        touch {output.fasta_clust_ref_file}
        """

rule analyse_blastn_raws:
    message:
        """
        analysing the results of the blastn false positive check for 
        {wildcards.sample} using R and taxonomizr package
        """
    input:
        blastfile = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_matches_nucleotide.m8",
        clusterfile = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads_clust.fasta.clstr"
    output:
        sampfinished = config["sub_dirs"]["raws_blastn_r"] + "/{sample}_finished.txt",
        viralreads = config["sub_dirs"]["raws_blastn_r"] + "/{sample}_readnames_virus_confirmed.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["raws_blastn_r"] + "/"
    log:
        "logs/" + config["sub_dirs"]["raws_blastn_r"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_blastn_r"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    resources:
        mem_mb=6000
    shell:
        """
        lengthblast=(`wc -l {input.blastfile}`) && \
        if [ ${{lengthblast}} -ge 2 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_raws_virus.R \
            --inputdiamond {input.blastfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.viralreads} \
            --programdir {params.basedir} --savdir {params.wrkdir} --clustfile {input.clusterfile}
        fi && \
        touch {output.sampfinished} && \
        touch {output.viralreads}
        """

rule summarise_raws:
    message:
        """
        analysing the results of the blastn false positive check for 
        {wildcards.sample} using R and taxonomizr package
        """
    input:
        viralreads = config["sub_dirs"]["raws_blastn_r"] + "/{sample}_readnames_virus_confirmed.txt",
        fastqreads = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads.fastq"
    output:
        fastqreadsfilt = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads_confirmed_blastn.fastq",
        sampfinished = config["sub_dirs"]["raws_results"] + "/{sample}_finished.txt",
        summarytable = config["sub_dirs"]["raws_results"] + "/{sample}Viral_hits_and_complexity.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["raws_results"] + "/",
        Viral_original_contigs = config["sub_dirs"]["Summary_results"] + "/{sample}_virus_hits_all_details_results_table.txt"
    log:
        "logs/" + config["sub_dirs"]["raws_results"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_results"] + "/{sample}.txt"
    conda: "R4_2_dada"
    threads: 1
    resources:
        mem_mb=3000
    shell:
        """
        lengthunassigned=(`wc -l {input.viralreads}`) && \
        if [ ${{lengthunassigned}} -ge 1 ]
        then
        cat {input.fastqreads} | grep --no-group-separator -A 3 -F -f "{input.viralreads}" > {output.fastqreadsfilt} && \
        Rscript {config[program_dir]}scripts/Summarise_raws.R \
            --inputvirusdetails {params.Viral_original_contigs} --name {params.samplename} \
            --threads {threads} --viral_readnames {input.viralreads} --fastqseqs {output.fastqreadsfilt} --output {output.summarytable} \
            --programdir {params.basedir} --savdir {params.wrkdir} && \
        touch {output.sampfinished}
        fi
        if [ ${{lengthunassigned}} -eq 0 ]
        then
        touch {output.sampfinished}
        touch {output.fastqreadsfilt}
        touch {output.summarytable}
        fi
        """

rule compile_raws_and_contigs:
    message:
        """
        Combining raw results with contig results, subsetting read names and preparing folders for individual species for 
        {wildcards.sample} using R
        """
    input:
        summarytable = config["sub_dirs"]["raws_results"] + "/{sample}Viral_hits_and_complexity.txt",
        fastqreads = config["sub_dirs"]["raws_blastn_check"] + "/{sample}_reads.fastq",
        contigsfile = config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa",
        summaryRdatafile = config["sub_dirs"]["Summary_results"] + "/{sample}_gather_summary_files_R_environment.Rdata"
    output:
        output_datatable = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_datatable.html",
        output_csv = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_sums.csv"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["compiled_summary"] + "/{sample}/",
        Viral_original_contigs = config["sub_dirs"]["Summary_results"] + "/{sample}_virus_hits_all_details_results_table.txt"
    log:
        "logs/" + config["sub_dirs"]["compiled_summary"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["compiled_summary"] + "/{sample}.txt"
    conda: "R4_2_dada"
    threads: 1
    resources:
        mem_mb=3000
    shell:
        """
        Rscript {config[program_dir]}scripts/Summarise_all_viral_reads_multi_methods.R \
            --inputvirusdetails {input.summarytable} --name {params.samplename} \
            --threads {threads} --contigs {input.contigsfile} --Rdatas {input.summaryRdatafile} \
            --programdir {params.basedir} --savdir {params.wrkdir} && \
        touch {output.output_datatable} && \
        touch {output.output_csv}
        """

rule extract_reads_and_contigs:
    message:
        """
        Extracting the final reads and contigs which aligned to each virus for 
        {wildcards.sample} .
        """
    input:
        R1 = host_removed_dataR1,
        R2 = host_removed_dataR2,
        output_datatable = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_datatable.html",
        contigfile = config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
    output:
        finished = config["sub_dirs"]["compiled_summary"] + "/{sample}/finished_extracting_reads.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["compiled_summary"] + "/{sample}"
    log:
        "logs/" + config["sub_dirs"]["compiled_summary"] + "/{sample}_extracting_reads.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["compiled_summary"] + "/{sample}_extracting_reads.txt"
    threads: 1
    resources:
        mem_mb=2000
    shell:
        """
        directories=$(ls -d {params.wrkdir}/*/) && \
        for dir in $directories; do
        if [[ ! ${{dir}} =~ _virusall_datatable_files ]]; then
        filtered_directories+=("${{dir}}")
        fi
        done
        for path in ${{filtered_directories[@]}}; do
        reads=$(ls ${{path}}*read_names*) && \
        if [ -n "${{reads}}" ]; then
        name=$(basename "${{reads}}" | sed 's/_read_names\.txt$//') && \
        zcat {input.R1} | grep -F -A 3 --no-group-separator -f ${{reads}} | gzip > ${{path}}${{name}}_reads_R1.fastq.gz && \
        zcat {input.R2} | grep -F -A 3 --no-group-separator -f ${{reads}} | gzip > ${{path}}${{name}}_reads_R2.fastq.gz && \
        contigs=$(find "${{path}}" -name "*contig_names*")
        fi && \
        if [ -n "${{contigs}}" ]; then
        cat ${{contigs}} | cut -f 1 > ${{contigs}}temp.txt && \
        seqkit grep -f ${{contigs}}temp.txt {input.contigfile} > ${{path}}${{name}}_contigs.fasta
        fi 
        done && \
        touch {output.finished}
        """