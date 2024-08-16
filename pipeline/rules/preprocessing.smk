# Write the end results (e.g., report or finalised clean reads. the end rule output)

rule record_start:
    output:
        "logs/start_time.txt"
    shell:
        """
        echo "Start time_"$(date) > {output}
        """

# wildcards.sample creates a wildcard variable for each of the returned values in the list of the 
# samples in the config file
# Can create the sample list automatically from just a directory and keyword as well
# Matt has sent an example from bacteria 

def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule filter_fastp:
    message:
        """
        preprocessing step 1 
        Trimming {wildcards.sample} for quality, adapters, and complexity using Fastp.
        """
    input:
        reads = get_raw_fastqs,
        date = "logs/start_time.txt"
    output:
        R1_filtered = config["sub_dirs"]["trim_dir"] + "/{sample}_R1.fastq.gz",
        R2_filtered = config["sub_dirs"]["trim_dir"] + "/{sample}_R2.fastq.gz",
        Readsfailed = config["sub_dirs"]["failed_trim_dir"] + "/{sample}_failed.fastq.gz",
        qualreporthtml = config["sub_dirs"]["trim_report_dir"] + "/{sample}_quality_report.html",
        qualreportjson = config["sub_dirs"]["trim_report_dir"] + "/{sample}_quality_report.json"
    params:
        minlen=config["minimum_length_filter_fastp"],
        complexity_threshold=config["complexity_threshold"],
        cut_front=config["front_window_cutsize"],
        minquality=config["min_qual_filter"],
        minavgqual=config["min_qual_filter_avg_read"],
        front_window_qual=config["min_qual_window"],
        adapters=config["Illumina_adapters"]
    log:
        "logs/" + config["sub_dirs"]["trim_dir"] + "/{sample}.log"
    threads: 8
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["trim_dir"] + "/{sample}.txt"
    shell:
        """
        fastp -i {input.reads[0]} -I {input.reads[1]} \
            -o {output.R1_filtered} -O {output.R2_filtered} \
            -j {output.qualreportjson} -h {output.qualreporthtml} \
            --adapter_fasta {params.adapters} \
            -q {params.minquality} \
            -e {params.minavgqual} \
            -w {threads} \
            -l {params.minlen} \
            --failed_out {output.Readsfailed} \
            --low_complexity_filter \
            --complexity_threshold {params.complexity_threshold} \
            --overrepresentation_analysis \
            -5 \
            -3 \
            --overrepresentation_sampling 100 \
            --cut_front_window_size {params.cut_front} \
            --cut_front_mean_quality {params.front_window_qual} \
            2> {log}
        """

rule filter_phiX:
    message:
        """
        preprocessing step 2 
        removing PhiX contamination from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["trim_dir"] + "/{sample}_R1.fastq.gz",
        R2 = config["sub_dirs"]["trim_dir"] + "/{sample}_R2.fastq.gz"
    output:
        R1 = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R1.fastq.gz",
        R2 = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R2.fastq.gz",
        samoutput = temp(config["sub_dirs"]["PhiX_dir"] + "/{sample}_PhiX.sam")
    params:
        Phixgenome = config["PhiX_genome_index"],
        unalignedreads = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R%.fastq.gz",
        alignedreads = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["PhiX_dir"] + "/{sample}.log"
    threads: 4
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["PhiX_dir"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.Phixgenome} -1 {input.R1} -2 {input.R2} \
            --un-conc-gz {params.unalignedreads} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        touch {output.R1} && \
        touch {output.R2}
        """

rule filter_CO1:
    message:
        """
        preprocessing step 3
        removing CO1 contamination from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R1.fastq.gz",
        R2 = config["sub_dirs"]["PhiX_dir"] + "/{sample}_R2.fastq.gz"
    output:
        unalignedreadsall1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        unalignedreadsall2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        unalignedsingles = config["sub_dirs"]["CO1_dir"] + "/{sample}_unaligned_singles.fastq.gz",
        alignedreads1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_aligned.fastq.gz",
        alignedreads2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_aligned.fastq.gz",
        samoutput = temp(config["sub_dirs"]["CO1_dir"] + "/{sample}_CO1.sam"),
        samoutputhits = temp(config["sub_dirs"]["CO1_dir"] + "/{sample}_CO1_hits.sam"),
        samoutputunaligned = temp(config["sub_dirs"]["CO1_dir"] + "/{sample}_CO1_unaligned_reads.sam")
    params:
        CO1genome = config["CO1_genome_index"],
        alignedreads = config["sub_dirs"]["CO1_dir"] + "/{sample}_R%_aligned.fastq.gz",
        samoutput = config["sub_dirs"]["CO1_dir"] + "/{sample}_CO1.sam"
    log:
        "logs/" + config["sub_dirs"]["CO1_dir"] + "/{sample}.log"
    threads: 4
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["CO1_dir"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.CO1genome} -1 {input.R1} -2 {input.R2} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        samtools view -@ {threads} -F 4 -h {output.samoutput} > {output.samoutputhits} && \
        samtools view -@ {threads} -f 4 -h {output.samoutput} > {output.samoutputunaligned} && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall1} -2 {output.unalignedreadsall2} -s {output.unalignedsingles} -n {output.samoutputunaligned} > /dev/null 2>&1
        """

rule filter_LSU:
    message:
        """
        preprocessing step 4 
        removing LSU contamination from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        R2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        unalignedsingles = config["sub_dirs"]["CO1_dir"] + "/{sample}_unaligned_singles.fastq.gz"
    output:
        unalignedsingles = config["sub_dirs"]["LSU_dir"] + "/{sample}_unaligned_singles.fastq.gz",
        unalignedreadsall1 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        unalignedreadsall2 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        alignedreads1 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        alignedreads2 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R2_aligned.fastq.gz",
        samoutput = temp(config["sub_dirs"]["LSU_dir"] + "/{sample}_LSU.sam"),
        samoutputhits = temp(config["sub_dirs"]["LSU_dir"] + "/{sample}_LSU_hits.sam"),
        samoutputunaligned = temp(config["sub_dirs"]["LSU_dir"] + "/{sample}_LSU_unaligned_reads.sam")
    params:
        LSUgenome = config["LSU_genome_index"],
        alignedreads = config["sub_dirs"]["LSU_dir"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["LSU_dir"] + "/{sample}.log"
    threads: 4
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["LSU_dir"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.LSUgenome} -1 {input.R1} -2 {input.R2} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        samtools view -@ {threads} -F 4 -h {output.samoutput} > {output.samoutputhits} && \
        samtools view -@ {threads} -f 4 -h {output.samoutput} > {output.samoutputunaligned} && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall1} -2 {output.unalignedreadsall2} -s {output.unalignedsingles} -n {output.samoutputunaligned} > /dev/null 2>&1
        """

rule filter_SSU:
    message:
        """
        preprocessing step 5 
        removing SSU contamination from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        R2 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        unalignedsingles = config["sub_dirs"]["LSU_dir"] + "/{sample}_unaligned_singles.fastq.gz"
    output:
        unalignedreadsall1 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        unalignedreadsall2 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        unalignedreadsall_singlesadded = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_unaligned_added_singles.fastq.gz",
        unalignedreadsall2_singlesadded = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_unaligned_added_singles.fastq.gz",
        alignedreads1 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        alignedreads2 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_aligned.fastq.gz",
        samoutput = temp(config["sub_dirs"]["SSU_dir"] + "/{sample}_SSU.sam"),
        samoutputhits = temp(config["sub_dirs"]["SSU_dir"] + "/{sample}_SSU_hits.sam"),
        unalignedsingles = config["sub_dirs"]["SSU_dir"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutputunaligned = temp(config["sub_dirs"]["SSU_dir"] + "/{sample}_SSU_unaligned_reads.sam")
    params:
        SSUgenome = config["SSU_genome_index"],
        alignedreads = config["sub_dirs"]["SSU_dir"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["SSU_dir"] + "/{sample}.log"
    threads: 4
    resources:
        mem_mb=12000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["SSU_dir"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.SSUgenome} -1 {input.R1} -2 {input.R2} -U {input.unalignedsingles} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        samtools view -@ {threads} -F 4 -h {output.samoutput} > {output.samoutputhits} && \
        samtools view -@ {threads} -f 4 -h {output.samoutput} > {output.samoutputunaligned} && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall1} -2 {output.unalignedreadsall2} -s {output.unalignedsingles} -n {output.samoutputunaligned}  > /dev/null 2>&1 && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall_singlesadded} -2 {output.unalignedreadsall2_singlesadded} -n {output.samoutputunaligned} > /dev/null 2>&1
        """
