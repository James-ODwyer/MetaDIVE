def LSU_hits(wildcards):
    if config["Rapid_host_id"] == "yes":
        return([
            config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt"
        ])
    elif config["Rapid_host_id"] == "no":
        return([
            config["sub_dirs"]["host_LCA_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt"
        ])

def SSU_hits(wildcards):
    if config["Rapid_host_id"] == "yes":
        return([
            config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt"
        ])
    elif config["Rapid_host_id"] == "no":
        return([
            config["sub_dirs"]["host_LCA_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt"
        ])

rule combine_barcode_marker_results_host:
    message:
        """
        host id step 1
        Extracting perfect(ish) CO1 allignments from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        CO1 = config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}_top_host_species_additional_stats_CO1.txt",
        LSU = LSU_hits,
        SSU = SSU_hits
    output:
        topspecies = config["sub_dirs"]["host_species_genomes"] + "/{sample}_top_host_species_overall.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_species_genomes"] + "/",
        CO1weight = config["CO1weight"],
        LSUweight = config["LSUweight"],
        SSUweight = config["SSUweight"]
    log:
        "logs/" + config["sub_dirs"]["host_species_genomes"] + "/id_{sample}.txt"
    threads: 1
    priority: 10
    conda: "Rdataplotting"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_species_genomes"] + "/id_{sample}.txt"
    shell:
        """
        Rscript {config[program_dir]}scripts/combine_barcode_host_ids.R \
            --inputCO1 {input.CO1} --inputLSU {input.LSU} --inputSSU {input.SSU} --name {params.samplename} \
            --programdir {params.basedir} --outdir {params.wrkdir} --CO1weighting {params.CO1weight} --LSUweighting {params.LSUweight} --SSUweighting {params.SSUweight} && \
        touch {output.topspecies}
        """

rule download_host_genome:
    message:
        """
        Downloading genome for host based on most probable host from the barcode markers CO1, LSU (rRNA subunit) and SSU (rRNA subunit)
        """
    input:
        topspecies = config["sub_dirs"]["host_species_genomes"] + "/{sample}_top_host_species_overall.txt"
    output:
        genome = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/top_host_species_CO1_genome.fa"
    params:
        programdir = config["program_dir"],
        outputfile = "top_host_species_CO1_genome.fa",
        outdir = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}",
        waittime = "20m",
        NCBI_key = config["NCBI_API_KEY"]
    log:
        "logs/" + config["sub_dirs"]["host_species_genomes"] + "/{sample}.txt"
    threads: 1
    priority: 10
    resources:
        genbank=1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_species_genomes"] + "/{sample}_download.txt"
    shell:
        """
        if [ ! -d {params.outdir} ]; then mkdir -p {params.outdir}; fi && \
        sbatch {config[program_dir]}scripts/downloadgenome.sh -i {input.topspecies} -o {params.outputfile} -n {params.NCBI_key} -p {params.outdir} -z {params.programdir} && \
        sleep {params.waittime} && \
        touch {output.genome}
        """

# May be worth changing the host index code so that it only looks for host species that haven't been assembled yet?

rule build_host_genome_idx:
    message:
        """
        Building genome index and aligning reads against the genome for removal
        """
    input:
        genome = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/top_host_species_CO1_genome.fa",
        R1 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_unaligned.fastq.gz",
        R2 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_unaligned.fastq.gz",
        unalignedsingles = config["sub_dirs"]["SSU_dir"] + "/{sample}_unaligned_singles.fastq.gz"
    output:
        unalignedreadsall1 = config["sub_dirs"]["host_remove_dir"] + "/{sample}_R1_host_removed.fastq.gz",
        unalignedreadsall2 = config["sub_dirs"]["host_remove_dir"] + "/{sample}_R2_host_removed.fastq.gz",
        unalignedreadsall_singlesadded = config["sub_dirs"]["host_remove_dir"] + "/{sample}_R1_host_removed_added_singles.fastq.gz",
        unalignedreadsall2_singlesadded = config["sub_dirs"]["host_remove_dir"] + "/{sample}_R2_host_removed_added_singles.fastq.gz",
        alignedreads1 = config["sub_dirs"]["host_remove_dir"] + "/aligned/{sample}_R1_aligned_host.fastq.gz",
        alignedreads2 = config["sub_dirs"]["host_remove_dir"] + "/aligned/{sample}_R2_aligned_host.fastq.gz",
        samoutput = temp(config["sub_dirs"]["host_remove_dir"] + "/aligned/{sample}_host.sam"),
        samoutputhits = temp(config["sub_dirs"]["host_remove_dir"] + "/aligned/{sample}_host_hits.sam"),
        unalignedsingles = config["sub_dirs"]["host_remove_dir"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutputunaligned= temp(config["sub_dirs"]["host_remove_dir"] + "/{sample}_host_unaligned_reads.sam")
    params:
        indexpath = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/index",
        host_genomeassembled = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/host_genome_complete.txt",
        host_genomeassembledgroup = config["sub_dirs"]["host_species_genomes"] + "/genomes/group/",
        host_genomeassembledgroupmarker = "/host_genome_complete.txt",
        temppathidx = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/index/Tempidx",
        unalignedreads = config["sub_dirs"]["host_remove_dir"] + "/{sample}_R%_host_removed.fastq.gz",
        alignedreads = config["sub_dirs"]["host_remove_dir"] + "/aligned/{sample}_R%_aligned_host.fastq.gz",
        nohostgenomemarker = "PROGRESS_STATUS/No_host_genome_was_found_{sample}.txt"
    log:
        "logs/" + config["sub_dirs"]["host_species_genomes"] + "/aligning/{sample}.txt"
    threads: 8
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_species_genomes"] + "/aligning/{sample}.txt"
    resources:
         mem_mb=10000
    shell:
        """
        if [ ! -d {params.indexpath} ]; then mkdir -p {params.indexpath}; fi && \
        length=(`wc -l {input.genome}`) && \
        if [ ${{length}} -ge 2 ]
        then
        bowtie2-build {input.genome} {params.temppathidx} --threads {threads} \
            --large-index -f
        bowtie2 -x {params.temppathidx} -1 {input.R1} -2 {input.R2} -U {input.unalignedsingles} \
            -p {threads} \
            --al-conc-gz {params.alignedreads} \
            -S {output.samoutput} \
            --fast-local\
            2> {log} && \
        samtools view -@ {threads} -F 4 -h {output.samoutput} > {output.samoutputhits} && \
        samtools view -@ {threads} -f 4 -h {output.samoutput} > {output.samoutputunaligned} && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall1} -2 {output.unalignedreadsall2} -s {output.unalignedsingles} -n {output.samoutputunaligned} > /dev/null 2>&1 && \
        samtools fastq -@ {threads} -1 {output.unalignedreadsall_singlesadded} -2 {output.unalignedreadsall2_singlesadded} -n {output.samoutputunaligned} > /dev/null 2>&1
        fi && \
        if [ ${{length}} -le 2 ]
        then
        touch {log}
        touch {output.samoutput}
        touch {output.samoutputhits}
        cp {input.R1} {output.unalignedreadsall1}
        cp {input.R2} {output.unalignedreadsall2}
        cp {input.unalignedsingles} {output.unalignedsingles}
        touch {output.alignedreads1}
        touch {output.alignedreads2}
        touch {output.unalignedsingles}
        touch {output.samoutputunaligned}
        touch {params.nohostgenomemarker}
        fi && \
        touch {output.alignedreads1} && \
        touch {output.alignedreads2} && \
        touch {output.alignedreads1} && \
        touch {output.alignedreads2} && \
        touch {output.unalignedsingles} && \
        touch {output.samoutputunaligned} && \
        touch {output.unalignedreadsall2_singlesadded} && \
        touch {output.unalignedreadsall_singlesadded} && \
        touch {log} && \
        touch {output.samoutput} && \
        touch {output.samoutputhits}
        """
#temp removed this part by inserting impossible length. Want to change to a blast index and blast assignment. It is far too prone to crashing as is due to 
# bowtie being so bad at large scaffolds

rule contigs_to_host:
    message:
        """
        Building genome index and aligning contigs from assembly against the host genome for removal
        """
    input:
        contigsfile = config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa",
        genome = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/top_host_species_CO1_genome.fa"
    output:
        unalignedcontigspre = temp(config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_removed_pre.fa"),
        unalignedcontigs = config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_removed.fa",
        alignedcontigs = config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_aligned_host.fa",
        samoutput = config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_alignment.sam",
        bamoutputhits = temp(config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_alignment_hits.bam"),
        samoutputnonhits = temp(config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_alignment_unaligned.sam")
    params:
        directory = config["sub_dirs"]["contig_dir_host_rem"],
        indexpath = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/index_blastn/",
        temppathidx = config["sub_dirs"]["host_species_genomes"] + "/genomes/{sample}/index_blastn/Tempidx",
        unalignedcontigsdir = config["sub_dirs"]["contig_dir_host_rem"],
        unalignedcontigsprename = config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs.fa"
    log:
        "logs/" + config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}.txt"
    threads: 4
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}.txt"
    resources:
         mem_mb=10000
    shell:
        """
        length=(`wc -l {input.genome}`) && \
        if [ ${{length}} -ge 1 ]; then
        minimap2 -ax asm5 --end-bonus 400 -s 450 --secondary=no -t {threads} {input.genome} {input.contigsfile} -o {output.samoutput} && \
        samtools view -f 4 -S {output.samoutput} -o {output.samoutputnonhits} && \
        samtools view -F 4 -b {output.samoutput} -o {output.bamoutputhits} && \
        reformat.sh -Xmx8g in={output.samoutputnonhits} out={output.unalignedcontigspre} && \
        seqkit rmdup -s < {output.unalignedcontigspre} > {output.unalignedcontigs}
        samtools fasta -@ {threads} {output.bamoutputhits} > {output.alignedcontigs}
        fi && \
        if [ ${{length}} -eq 0 ]; then
        cp {input.contigsfile} {params.unalignedcontigsdir} && \
        mv {params.unalignedcontigsprename} {output.unalignedcontigs} && \
        touch {output.alignedcontigs} && \
        touch {output.samoutput}
        fi && \
        touch {output.alignedcontigs} && \
        touch {output.unalignedcontigs} && \
        touch {output.unalignedcontigspre} && \
        touch {output.samoutput} && \
        touch {output.bamoutputhits} && \
        touch {output.samoutputnonhits}
        """

rule align_raws_to_host_contigs:
    message:
        """
        Assembly step 4 
        Aligning original reads to generated contigs for {wildcards.sample} using Bowtie2
        """
    input:
        alignedcontigs = config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_aligned_host.fa",
        reads1 = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R1.fastq.gz",
        reads2 = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_R2.fastq.gz"
    output:
        unalignedreads1 = config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R1.fastq.gz",
        unalignedreads2 = config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R2.fastq.gz",
        alignedreads1 = temp(config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R1_aligned.fastq.gz"),
        alignedreads2 = temp(config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R2_aligned.fastq.gz"),
        samoutput = temp(config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_SSU.sam"),
        samoutputhits = config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_host_contig_hits.sam"
    params:
        indexpath = config["sub_dirs"]["temp_genome_path"] + "/host_contigs/{sample}/",
        temppathidx = config["sub_dirs"]["temp_genome_path"] + "/host_contigs/{sample}/{sample}_Tempidx",
        unalignedreads = config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R%.fastq.gz",
        alignedreads = config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_R%_aligned.fastq.gz"
    log:
        "logs/" + config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}.log"
    threads: 8
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}.txt"
    shell:
        """
        length=(`wc -l {input.alignedcontigs}`) && \
        if [ ${{length}} -ge 2 ]
        then
        if [ ! -d {params.indexpath} ]; then mkdir -p {params.indexpath}; fi && \
        bowtie2-build {input.alignedcontigs} {params.temppathidx} --threads {threads} \
            --large-index
        bowtie2 -x {params.temppathidx} -1 {input.reads1} -2 {input.reads2} \
            --un-conc-gz {params.unalignedreads} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --fast \
            2> {log} && \
        samtools view -F 4 {output.samoutput} > {output.samoutputhits}
        fi && \
        if [ ${{length}} -le 2 ]
        then
        cp {input.reads1} {output.unalignedreads1}
        cp {input.reads2} {output.unalignedreads2}
        fi && \
        touch {output.unalignedreads1} && \
        touch {output.unalignedreads2} && \
        touch {output.alignedreads1} && \
        touch {output.alignedreads2} && \
        touch {output.samoutput} && \
        touch {output.samoutputhits}
        """
