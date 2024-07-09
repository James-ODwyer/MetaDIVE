# want to 1. Extract all Viral genome hits per dramv annotation
# 2. link and download docsum for relevant virus genomes
# 3. make sure each genome is only present once
# 4. Download genomes
# 5. align raw reads to virus genome and call % coverage
# 6. generate summary of what % of each genome each virus was found in

def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule find_viral_genomes:
    message:
        """
        Advanced viral Identification step 5.  
        Downloading the genomes of all viruses identified as references in the dram annotation
        """
    input:
        finishedcontigassignment_pipeline = config["sub_dirs"]["finished"] + "/{sample}_finished"
    output:
        finished = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}/finished.txt"
    params:
        programdir = config["program_dir"],
        topviralcontigs = config["sub_dirs"]["Summary_results"] + "/{sample}_top10_Taxids_Viral_contigs.txt",
        outdir = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}",
        NCBI_key = config["NCBI_API_KEY"]
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present"] + "/{sample}.txt"
    threads: 1
    resources:
         genbank=1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present"] + "/{sample}_download.txt"
    shell:
        """
        if [ ! -d {params.outdir} ]; then mkdir -p {params.outdir}; fi && \
        length2=(`wc -l {params.topviralcontigs}`) && \
        if [[ $length2 -ge 1 ]]; then
        sbatch {config[program_dir]}scripts/download_viral_genomes.sh -j {params.topviralcontigs} -n {params.NCBI_key} -p {params.outdir} -z {params.programdir} && \
        while [ ! -e "{output.finished}" ]
        do
        sleep 180 # wait for 3 min before checking again
        done
        fi && \
        touch {output.finished}
        """

rule generate_viral_genome_indexes:
    message:
        """
        Advanced viral Identification step 6.  
        Creating bowtie indexes of reference genomes and aligning raw reads to each genome
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}/finished.txt",
        unalignedreads1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        unalignedreads2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_unaligned.fastq.gz"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}/finished_align.txt"
    params:
        programdir = config["program_dir"],
        samplename = "{sample}",
        ingenomedir = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}",
        inhostremoveddir = config["sub_dirs"]["CO1_dir"],
        outconsensusdir = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}"
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present_summary"] + "/{sample}.txt"
    threads: 4
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present"] + "/{sample}_assemble.txt"
    shell:
        """
        if [ ! -d {params.outconsensusdir} ]; then mkdir -p {params.outconsensusdir}; fi && \
        length=(`ls {params.ingenomedir}/* | wc -l`) && \
        if [ ${{length}} -ge 2 ]
        then
        bash {config[program_dir]}scripts/assemble_viral_genomes.sh -g {params.ingenomedir} -j {input.unalignedreads1} -k {input.unalignedreads2} -p {params.inhostremoveddir} -r {params.outconsensusdir} -t {params.samplename} -z {params.programdir} && \
        while [ ! -e "{output.finished2}" ]
        do
        sleep 180 # wait for 3 min before checking again
        done
        fi && \
        if [ ${{length}} -le 1 ]
        then
        echo " No genes were detected for sample {params.samplename} " && \ 
        touch {output.finished2} 
        fi && \
        touch {output.finished2}
        """

rule filter_assembled_viral_genomes:
    message:
        """
        Advanced viral Identification step 7.  
        Filtering genomes generated down to only include genomes of 50% or higher coverage
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}/finished_align.txt"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present_summary"] + "/{sample}/finished2.txt"
    params:
        indir = config["sub_dirs"]["Viral_genomes_present"] + "/{sample}",
        outdir = config["sub_dirs"]["Viral_genomes_present_summary"] + "/{sample}",
        samp = "{sample}"
    threads: 1
    shell:
        """
        if [ ! -d {params.outdir} ] 
        then 
        mkdir -p {params.outdir} 
        fi && \
        bash {config[program_dir]}scripts/subset_viral_genomes_50cov.sh -i {params.indir} -o {output.finished2} -p {params.outdir} -z {params.samp} && \
        while [ ! -e "{output.finished2}" ]
        do
        sleep 180 # wait for 3 min before checking again
        done
        touch {output.finished2}
        """

rule generate_viral_genomes_guided_alignments:
    message:
        """
        Advanced viral Identification step 6.  
        Creating bowtie indexes of reference genomes and aligning raw reads to each genome
        """
    input:
        finished = config["sub_dirs"]["Viral_genomes_present_summary"] + "/{sample}/finished2.txt",
        unalignedreads1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        unalignedreads2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_unaligned.fastq.gz"
    output:
        finished2 = config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}/finished2.txt"
    params:
        programdir = config["program_dir"],
        samplename = "{sample}",
        ingenomedir = config["sub_dirs"]["Viral_genomes_present_summary"] + "/{sample}",
        inhostremoveddir = config["sub_dirs"]["SSU_dir"],
        outconsensusdir = config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}",
        NCBI_key = config["NCBI_API_KEY"]
    log:
        "logs/" + config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}.txt"
    threads: 1
    resources:
         genbank=1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}_download.txt"
    shell:
        """
        if [ ! -d {params.outconsensusdir} ]; then mkdir -p {params.outconsensusdir}; fi && \
        length=(`ls {params.ingenomedir}/* | wc -l`) && \
        if [ ${{length}} -ge 3 ]
        then
        sbatch {config[program_dir]}scripts/download_variety_genomes_high_completion_viruses.sh -g {params.ingenomedir} -j {input.unalignedreads1} -k {input.unalignedreads2} -n {params.NCBI_key} -p {params.inhostremoveddir} -r {params.outconsensusdir} -t {params.samplename} -z {params.programdir} && \
        while [ ! -e "{output.finished2}" ]
        do
        sleep 60 # wait for 3 min before checking again
        done
        fi && \
        if [ ${{length}} -le 2 ]
        then
        echo " No genes were detected for sample {params.samplename} " && \ 
        touch {output.finished2} 
        fi && \
        touch {output.finished2}
        """
