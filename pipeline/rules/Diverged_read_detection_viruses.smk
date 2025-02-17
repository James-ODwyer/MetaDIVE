
def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]


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
            "30"
        ])
    if config["Sensitivity"] == "Medium":
        return([
            "24"
        ])
    if config["Sensitivity"] == "Medium-Low":
        return([
            "18"
        ])
    if config["Sensitivity"] == "Low":
        return([
            "16"
        ])
    if config["Sensitivity"] == "Ultra_low":
        return([
            "16"
        ])


def Diamond_memrawMB(wildcards):
    if config["Sensitivity"] == "VHigh":
        return "36000"
    if config["Sensitivity"] == "High":
        return "30000"
    if config["Sensitivity"] == "Medium":
        return "24000"
    if config["Sensitivity"] == "Medium-Low":
        return "18000"
    if config["Sensitivity"] == "Low":
        return "16000"
    if config["Sensitivity"] == "Ultra_low":
        return "16000"

rule prep_contigs_plus_reads_rediamond:
    message:
        """
        Assembly step 1 divergent reads detection 
        removing identified reads and contigs from remaining data for {wildcards.sample}
        """
    input:
        unalignedreads1 = raws_after_hostR1,
        unalignedreads2 = raws_after_hostR2,
        unalignedcontigs = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa",
        identified_raws = config["sub_dirs"]["raws_blastn_r"] + "/{sample}_readnames_virus_confirmed.txt",
        output_datatable = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_datatable.html"
    output:
        fastqreadsfiltR1 = temp(config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_R1_virdepl.fastq"),
        fastqreadsfiltR2 = temp(config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_R2_virdepl.fastq"),
        interleavedfiles = temp(config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_interleaved.fastq"),
        interleavedfasta = temp(config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_interleaved.fasta"),
        interleavedfiles_both = temp(config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_interleaved_both.fasta"),
        diamondfile = config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_matches.m8"
    params:
        databasenr = config["diamond_database"],
        unalignedreads1_prior = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1.fastq.gz",
        unalignedreads2_prior = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2.fastq.gz",
        taxidlist = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_detected_viral_taxids.csv",
        diamondsensitivity = config["Divergent_reads_and_contigs_sensitivity"],
        diamondmem = Diamond_memraw
    log:
        "logs/" + config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}.log"
    threads: 4
    priority: 5
    conda: "snakemake7"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}.txt"
    resources:
        mem_mb = Diamond_memrawMB
    shell:
        """
        seqkit grep -v -f "{input.identified_raws}" "{input.unalignedreads1}" > "{output.fastqreadsfiltR1}" && \
        seqkit grep -v -f "{input.identified_raws}" "{input.unalignedreads2}" > "{output.fastqreadsfiltR2}" && \
        seqfu ilv -1 {output.fastqreadsfiltR1} -2 {output.fastqreadsfiltR2} -o {output.interleavedfiles} && \
        seqtk seq -A {output.interleavedfiles} > {output.interleavedfasta} && \
        cat {output.interleavedfasta} {input.unalignedcontigs} > {output.interleavedfiles_both} && \
        diamond blastx --db {params.databasenr} \
            --query {output.interleavedfiles} \
             --iterate \
            --{params.diamondsensitivity} \
            --taxonlist "$(paste -sd, {params.taxidlist})" \
            --max-target-seqs 1 \
            -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp \
            --evalue 0.001 \
            --threads {threads} \
            -o {output.diamondfile} \
            --memory-limit {params.diamondmem} \
            2> {log}
        """

rule analyse_blastx_sensitive:
    message:
        """
        analysing the results of the blastn false positive check for 
        {wildcards.sample} using R and taxonomizr package
        """
    input:
        diamondfile = config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_matches.m8",
        resultscsv_virus_all = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_sums.csv"
    output:
        sampfinished = config["sub_dirs"]["compiled_summary"] + "/{sample}_diverged_finished.txt",
        resultscsv_virus_incl_diverged = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_sums_incl_diverged_reads.csv"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["compiled_summary"] + "/{sample}/"
    log:
        "logs/" + config["sub_dirs"]["raws_blastn_r"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_blastn_r"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 2
    resources:
        mem_mb=4000
    shell:
        """
        lengthblast=(`wc -l {input.diamondfile}`) && \
        if [ ${{lengthblast}} -ge 2 ]
        then
        Rscript {config[program_dir]}scripts/add_diverged_reads_and_contigs_counts.R \
            --inputdiamond {input.diamondfile} --name {params.samplename} --inputvirusall {input.resultscsv_virus_all} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.resultscsv_virus_incl_diverged} \
            --programdir {params.basedir} --savdir {params.wrkdir}
        fi && \
        touch {output.sampfinished} && \
        touch {output.resultscsv_virus_incl_diverged}
        """

rule extract_diverged_reads_and_contigs:
    message:
        """
        Extracting the final reads and contigs which aligned to each virus for 
        {wildcards.sample} .
        """
    input:
        fastqreadsfiltR1 = config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_R1_virdepl.fastq",
        fastqreadsfiltR2 = config["sub_dirs"]["Viral_diamond_sensitive_check"] + "/{sample}_R2_virdepl.fastq",
        unalignedcontigs = config["sub_dirs"]["contigs_assigned"] + "/{sample}_matches_unassigned.fa",
        resultscsv_virus_incl_diverged = config["sub_dirs"]["compiled_summary"] + "/{sample}/{sample}_virusall_sums_incl_diverged_reads.csv",
        contigfile = config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
    output:
        finished = config["sub_dirs"]["compiled_summary"] + "/{sample}/finished_extracting_reads_diverged.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["compiled_summary"] + "/{sample}"
    log:
        "logs/" + config["sub_dirs"]["compiled_summary"] + "/{sample}_extracting_reads_diverged.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["compiled_summary"] + "/{sample}_extracting_reads_diverged.txt"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        """
        filtered_directories=()
        for dir in {params.wrkdir}/*/; do
            [[ ${{dir}} != *_virusall_datatable_files/ ]] && filtered_directories+=("${{dir}}")
        done && \
        for path in "${{filtered_directories[@]}}"; do
            contigs=("${{path}}"*divergent*)
            if [[ -f "${{contigs[0]}}" ]]; then
                read_file="${{contigs[0]}}"
                cut -f 1 "${{contigs}}" > "${{contigs}}temp.txt"
                name=${{read_file%_divergent_reads_and_contigs.txt}}
                name=${{name##*/}}
                cat {input.fastqreadsfiltR1} | grep -F -A 3 --no-group-separator -f "${{contigs}}temp.txt" | gzip > "${{path}}${{name}}_diverged_reads_R1.fastq.gz"
                cat {input.fastqreadsfiltR2} | grep -F -A 3 --no-group-separator -f "${{contigs}}temp.txt" | gzip > "${{path}}${{name}}_diverged_reads_R2.fastq.gz"
                seqkit grep -f "${{contigs}}temp.txt" {input.unalignedcontigs} > "${{path}}${{name}}_diverged_contigs.fasta"
            fi
        done && \
        touch {output.finished}
        """
