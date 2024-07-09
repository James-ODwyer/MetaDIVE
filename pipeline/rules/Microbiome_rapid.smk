rule strong_filter_LSU_rapid:
    message:
        """
        host id step 2
        Extracting perfect(ish) LSU allignments from {wildcards.sample} using bowtie2 local alignments for rapid analysis
        """
    input:
        R1 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        R2 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R2_aligned.fastq.gz"
    output:
        R1poshitssub = temp(config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_R1_subset.fq.gz"),
        R2poshitssub = temp(config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_R2_subset.fq.gz"),        
        alignedreads1 = temp(config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_R1_aligned.fastq.gz"),
        alignedreads2 = temp(config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_R2_aligned.fastq.gz"),
        samoutput = temp(config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_LSU.sam"),
        sampos = config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_pos_hits.sam"
    params:
        LSUgenome = config["LSU_genome_index"],
       #unaligned would go here but no point in having it for this step of the analysis
        alignedreads = config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_R%_aligned.fastq.gz",
        # Make samoutput a temp file
    log:
        "logs/" + config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}.txt"
    threads: 8
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}.txt"
    shell:
        """
        seqtk sample -s2468 {input.R1} 50000 > {output.R1poshitssub} && \
        seqtk sample -s2468 {input.R2} 50000 > {output.R2poshitssub} && \
        bowtie2 -x {params.LSUgenome} -1 {output.R1poshitssub} -2 {output.R2poshitssub} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --sensitive \
            2> {log} && \
        samtools view -F 4 -h {output.samoutput} > {output.sampos}
        """

rule strong_filter_SSU_rapid:
    message:
        """
        host id step 3
        Extracting perfect(ish) SSU allignments from {wildcards.sample} using bowtie2 local alignments for rapid analysis
        """
    input:
        R1 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        R2 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_aligned.fastq.gz"
    output:
        R1poshitssub = temp(config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_R1_subset.fq.gz"),
        R2poshitssub = temp(config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_R2_subset.fq.gz"),        
        alignedreads1 = temp(config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_R1_aligned.fastq.gz"),
        alignedreads2 = temp(config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_R2_aligned.fastq.gz"),
        sampos = config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_pos_hits.sam",
        samoutput = temp(config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_SSU.sam")
    params:
        SSUgenome = config["SSU_genome_index"],
       #unaligned would go here but no point in having it for this step of the analysis
        alignedreads = config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_R%_aligned.fastq.gz",
        # Make samoutput a temp file
    log:
        "logs/" + config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}.txt"
    threads: 8
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}.txt"
    shell:
        """
        seqtk sample -s2468 {input.R1} 50000 > {output.R1poshitssub}
        seqtk sample -s2468 {input.R2} 50000 > {output.R2poshitssub}
        bowtie2 -x {params.SSUgenome} -1 {output.R1poshitssub} -2 {output.R2poshitssub} \
            --al-conc-gz {params.alignedreads} \
            -p {threads} \
            -S {output.samoutput} \
            --sensitive \
            2> {log} && \
        samtools view -F 4 -h {output.samoutput} > {output.sampos}
        """

rule align_all_returned_LSU_rapid:
    message:
        """
        microbiome step 2
        Using mmseq2 with LCA based methods to classify taxonomy for rapid analysis
        """
    input:
        poshitssam = config["sub_dirs"]["host_detect_LSU_rapid"] + "/{sample}_pos_hits.sam"
    output:
        R1poshits = temp(config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_R1_pos_hits.fq"),
        R2poshits = temp(config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_R2_pos_hits.fq"),
        singletons = temp(config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_singletons_pos_hits.fq"),
        taxresulttsv = config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_taxonomy_report_LSU_rapid.tsv",
        taxresultkraken = config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_taxonomy_report_LSU_kraken_rapid.tsv"
    params:
        LSUgenome = config["LSU_genome_index_mmseq"],
        tempqueryDBfile = config["program_dir"] + "mmeseq/TEMP/LSU/rapid/{sample}/{sample}_TEMP_DB",
        tempqueryDBlocation = config["program_dir"] + "mmeseq/TEMP/LSU/rapid/{sample}",
        taxresult = config["program_dir"] + "mmeseq/TEMP/LSU/rapid/{sample}/{sample}_taxonomy_report"
    log:
        "logs/" + config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}.log"
    threads: 16
    priority: 10
    resources:
         mem_mb=60000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.tempqueryDBlocation} ] 
        then 
        mkdir -p {params.tempqueryDBlocation} 
        fi && \
        length=(`wc -l {input.poshitssam}`) && \
        if [ ${{length}} -ge 1 ]
        then
        samtools fastq -n {input.poshitssam} -1 {output.R1poshits} -2 {output.R2poshits} -s {output.singletons} && \
        mmseqs createdb {output.R1poshits} {output.R2poshits} {params.tempqueryDBfile} --dbtype 2 && \
        mmseqs taxonomy {params.tempqueryDBfile} {params.LSUgenome} {params.taxresult} {params.tempqueryDBlocation} \
            -s 1.0 --max-seqs 50 --threads {threads} --search-type 3 --tax-lineage 1 --split-memory-limit 22G --remove-tmp-files 1 && \
        mmseqs createtsv {params.tempqueryDBfile} {params.taxresult} {output.taxresulttsv} && \
        mmseqs taxonomyreport {params.LSUgenome} {params.taxresult} {output.taxresultkraken} --threads {threads} && \
        rm -r {params.tempqueryDBlocation}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.R1poshits} && \
        touch {output.R2poshits} && \
        touch {output.singletons} && \
        touch {output.taxresulttsv} && \
        touch {output.taxresultkraken}
        """

rule align_all_returned_SSU_rapid:
    message:
        """
        microbiome step 2
        Using mmseq2 with LCA based methods to classify taxonomy for a rapid subset of the SSU data (host id) for rapid analysis
        """
    input:
        poshitssam = config["sub_dirs"]["host_detect_SSU_rapid"] + "/{sample}_pos_hits.sam"
    output:
        R1poshits = temp(config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_R1_pos_hits.fq"),
        R2poshits = temp(config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_R2_pos_hits.fq"),
        singletons = temp(config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_singletons_pos_hits.fq"),
        taxresulttsv = config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_taxonomy_report_SSU_rapid.tsv",
        taxresultkraken = config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_taxonomy_report_SSU_kraken_rapid.tsv"
    params:
        SSUgenome = config["SSU_genome_index_mmseq"],
        tempqueryDBfile = config["program_dir"] + "mmeseq/TEMP/SSU/rapid/{sample}/{sample}_TEMP_DB",
        tempqueryDBlocation = config["program_dir"] + "mmeseq/TEMP/SSU/rapid/{sample}",
        taxresult = config["program_dir"] + "mmeseq/TEMP/SSU/rapid/{sample}/{sample}_taxonomy_report"
    log:
        "logs/" + config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}.log"
    threads: 16
    priority: 10
    resources:
         mem_mb=60000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.tempqueryDBlocation} ] 
        then 
        mkdir -p {params.tempqueryDBlocation} 
        fi && \
        length=(`wc -l {input.poshitssam}`) && \
        if [ ${{length}} -ge 1 ]
        then
        samtools fastq -n {input.poshitssam} -1 {output.R1poshits} -2 {output.R2poshits} -s {output.singletons} && \
        mmseqs createdb {output.R1poshits} {output.R2poshits} {params.tempqueryDBfile} --dbtype 2 && \
        mmseqs taxonomy {params.tempqueryDBfile} {params.SSUgenome} {params.taxresult} {params.tempqueryDBlocation} \
            -s 1.0 --max-seqs 50 --threads {threads} --search-type 3 --tax-lineage 1 --split-memory-limit 22G --remove-tmp-files 1 && \
        mmseqs createtsv {params.tempqueryDBfile} {params.taxresult} {output.taxresulttsv} && \
        mmseqs taxonomyreport {params.SSUgenome} {params.taxresult} {output.taxresultkraken} --threads {threads} && \
        rm -r {params.tempqueryDBlocation}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.R1poshits} && \
        touch {output.R2poshits} && \
        touch {output.singletons} && \
        touch {output.taxresulttsv} && \
        touch {output.taxresultkraken}
        """

rule Generate_figures_LSU_mapping_rapid:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from LSU rapid subset
        """
    input:
        taxresultkraken = config["sub_dirs"]["host_assignLCA_rapid_LSU"] + "/{sample}_taxonomy_report_LSU_kraken_rapid.tsv"
    output:
        sampfinished = config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/",
        nHost = config["Hostdetect"],
        marker = "LSU"
    log:
        "logs/" + config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    priority: 10
    shell:
        """
        Rscript {config[program_dir]}scripts/Generate_LCAs_barcode_markers.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} && \
        touch {output.sampfinished}
        """

rule Generate_figures_SSU_mapping_rapid:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from SSU rapid subset 
        """
    input:
        taxresultkraken = config["sub_dirs"]["host_assignLCA_rapid_SSU"] + "/{sample}_taxonomy_report_SSU_kraken_rapid.tsv"
    output:
        sampfinished = config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/",
        nHost = config["Hostdetect"],
        marker = "SSU"
    log:
        "logs/" + config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    priority: 10
    shell:
        """
        Rscript {config[program_dir]}scripts/Generate_LCAs_barcode_markers.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} && \
        touch {output.sampfinished}
        """

rule collate_barcodes_results_rapid:
    message:
        """
        microbiome step 4
        collate results from each barcod marker rapid analysis
        """
    input:
        sampfinishedSSUrapid = config["sub_dirs"]["host_LCA_rapid_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt",
        sampfinishedLSUrapid = config["sub_dirs"]["host_LCA_rapid_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt",
        sampfinishedCO1 = config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}_top_host_species_additional_stats_CO1.txt"
    output:
        analysisfinished =  config["sub_dirs"]["progress_barcodes"] + "/rapid/{sample}_completed_fullbarcoding_analysis.txt"
    threads: 1
    priority: 10
    shell:
        """
        touch {output.analysisfinished}
        """