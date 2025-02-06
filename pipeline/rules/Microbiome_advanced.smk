
rule strong_filter_CO1:
    message:
        """
        host id step 1
        Extracting perfect(ish) CO1 allignments from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R1_aligned.fastq.gz",
        R2 = config["sub_dirs"]["CO1_dir"] + "/{sample}_R2_aligned.fastq.gz"
    output:
        R1spos = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R2s.fastq",
        unalignedsingles = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutput = temp(config["sub_dirs"]["host_detect_CO1"] + "/{sample}_CO1.sam"),
        sampos = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_pos_hits.sam",
        samposnohead= config["sub_dirs"]["host_detect_CO1"] + "/{sample}_pos_hits_no_head.sam"
    params:
        CO1genome = config["CO1_genome_index"]
    log:
        "logs/03_CO1/HOSTREMOVAL/01_STRONG_ALIGN/{sample}.log"
    threads: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_detect_CO1"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.CO1genome} -1 {input.R1} -2 {input.R2} \
            -p {threads} \
            -S {output.samoutput} \
            --very-fast \
            2> {log} && \
        samtools view -@ {threads} -h -F 4 {output.samoutput} > {output.sampos} && \
        samtools view -@ {threads} -F 4 {output.samoutput} > {output.samposnohead} && \
        samtools fastq -1 {output.R1spos} -2 {output.R2spos} -s {output.unalignedsingles} -n {output.sampos} > /dev/null 2>&1
        """

rule generate_CO1_contigs:
    message:
        """
        Assembly step 1 
        Generating contigs from {wildcards.sample} using Megahit kmer based de novo assembly
        """
    input:
        R1s = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R1s.fastq",
        R2s = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R2s.fastq"
    output:
        contigsfile = config["sub_dirs"]["host_contigs_CO1"] + "/{sample}/final.contigs.fa"
    params:
        megahitpath = config["sub_dirs"]["host_contigs_CO1"] +"/{sample}"
    log:
        "logs/" + config["sub_dirs"]["host_contigs_CO1"] + "/{sample}.log"
    threads: 4
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_contigs_CO1"] + "/{sample}.txt"
    resources:
        mem_mb=6000
    shell:
        """
        megahit -1 {input.R1s} -2 {input.R2s} \
            --presets meta-sensitive \
            --memory 0.95 --force \
            -t {threads} \
            -o {params.megahitpath} \
            --min-contig-len 600 \
            2> {log}
        """

rule strong_filter_LSU:
    message:
        """
        host id step 2
        Extracting perfect(ish) LSU allignments from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        R2 = config["sub_dirs"]["LSU_dir"] + "/{sample}_R2_aligned.fastq.gz",
        samoutputCO1 = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_pos_hits.sam"
    output:
        R1spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R2s.fastq",
        unalignedsingles = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutput = temp(config["sub_dirs"]["host_detect_LSU"] + "/{sample}_LSU.sam"),
        sampos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_pos_hits.sam",
        samposnohead = temp(config["sub_dirs"]["host_detect_LSU"] + "/{sample}_pos_hits_no_head.sam")
    params:
        LSUgenome = config["LSU_genome_index"]
       #unaligned would go here but no point in having it for this step of the analysis
    log:
        "logs/04_LSU/HOSTREMOVAL/01_STRONG_ALIGN/{sample}.log"
    threads: 2
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_detect_LSU"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.LSUgenome} -1 {input.R1} -2 {input.R2} \
            -p {threads} \
            -S {output.samoutput} \
            --very-fast \
            2> {log} && \
        samtools view -@ {threads} -h -F 4 {output.samoutput} > {output.sampos} && \
        samtools view -@ {threads} -F 4 {output.samoutput} > {output.samposnohead} && \
        samtools fastq -1 {output.R1spos} -2 {output.R2spos} -s {output.unalignedsingles} -n {output.sampos} > /dev/null 2>&1
        """

rule generate_LSU_contigs:
    message:
        """
        Assembly step 1 
        Generating contigs from {wildcards.sample} using Megahit kmer based de novo assembly
        """
    input:
        R1spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R2s.fastq"
    output:
        contigsfile = config["sub_dirs"]["host_contigs_LSU"] + "/{sample}/final.contigs.fa"
    params:
        megahitpath = config["sub_dirs"]["host_contigs_LSU"] +"/{sample}"
    log:
        "logs/" + config["sub_dirs"]["host_contigs_LSU"] + "/{sample}.log"
    threads: 4
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_contigs_LSU"] + "/{sample}.txt"
    resources:
        mem_mb=6000
    shell:
        """
        megahit -1 {input.R1spos} -2 {input.R2spos} \
            --presets meta-sensitive \
            --memory 0.95 --force \
            -t {threads} \
            -o {params.megahitpath} \
            --min-contig-len 600 \
            2> {log}
        """

rule strong_filter_SSU:
    message:
        """
        host id step 3
        Extracting perfect(ish) SSU allignments from {wildcards.sample} using bowtie2 local alignments
        """
    input:
        R1 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R1_aligned.fastq.gz",
        R2 = config["sub_dirs"]["SSU_dir"] + "/{sample}_R2_aligned.fastq.gz",
        samposLSU = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_pos_hits.sam"
    output:
        R1spos = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R2s.fastq",
        unalignedsingles = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_unaligned_singles.fastq.gz",
        samoutput = temp(config["sub_dirs"]["host_detect_SSU"] + "/{sample}_SSU.sam"),
        samposnohead = temp(config["sub_dirs"]["host_detect_SSU"] + "/{sample}_pos_hits_no_head.sam"),
        sampos = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_pos_hits.sam"
    params:
        SSUgenome = config["SSU_genome_index"]
    log:
        "logs/05_SSU/HOSTREMOVAL/01_STRONG_ALIGN/{sample}.log"
    threads: 2
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_detect_SSU"] + "/{sample}.txt"
    shell:
        """
        bowtie2 -x {params.SSUgenome} -1 {input.R1} -2 {input.R2} \
            -p {threads} \
            -S {output.samoutput} \
            --very-fast \
            2> {log} && \
        samtools view -@ {threads} -h -F 4 {output.samoutput} > {output.sampos} && \
        samtools view -@ {threads} -F 4 {output.samoutput} > {output.samposnohead} && \
        samtools fastq -1 {output.R1spos} -2 {output.R2spos} -s {output.unalignedsingles} -n {output.sampos} > /dev/null 2>&1
        """

rule generate_SSU_contigs:
    message:
        """
        Assembly step 1 
        Generating contigs from {wildcards.sample} using Megahit kmer based de novo assembly
        """
    input:
        R1s = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R1s.fastq",
        R2s = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R2s.fastq"
    output:
        contigsfile = config["sub_dirs"]["host_contigs_SSU"] + "/{sample}/final.contigs.fa"
    params:
        megahitpath = config["sub_dirs"]["host_contigs_SSU"] +"/{sample}"
    log:
        "logs/" + config["sub_dirs"]["host_contigs_SSU"] + "/{sample}.log"
    threads: 4
    priority: 10
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_contigs_SSU"] + "/{sample}.txt"
    resources:
        mem_mb=6000
    shell:
        """
        megahit -1 {input.R1s} -2 {input.R2s} \
            --presets meta-sensitive \
            --memory 0.95 --force \
            -t {threads} \
            -o {params.megahitpath} \
            --min-contig-len 600 \
            2> {log}
        """

rule align_all_returned_CO1:
    message:
        """
        microbiome step 2
        Using mmseq2 with LCA based methods to classify taxonomy
        """
    input:
        R1s = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R1s.fastq",
        R2s = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_R2s.fastq",
        samposnohead = config["sub_dirs"]["host_detect_CO1"] + "/{sample}_pos_hits_no_head.sam"
    output:
        R1poshitssub = temp(config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_R1_pos_hits_subset.fq"),
        R2poshitssub = temp(config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_R2_pos_hits_subset.fq"),
        taxresulttsv = config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_taxonomy_report_CO1.tsv",
        taxresultkraken = config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_taxonomy_report_CO1_kraken.tsv"
    params:
        CO1genome = config["CO1_genome_index_mmseq"],
        tempqueryDBfile = config["program_dir"] + "mmeseq/TEMP/CO1/{sample}/{sample}_TEMP_DB",
        tempqueryDBlocation = config["program_dir"] + "mmeseq/TEMP/CO1/{sample}",
        taxresult = config["program_dir"] + "mmeseq/TEMP/CO1/{sample}/{sample}_taxonomy_report"
    log:
        "logs/" + config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}.log"
    threads: 4
    resources:
         mem_mb=24000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.tempqueryDBlocation} ] 
        then 
        mkdir -p {params.tempqueryDBlocation} 
        fi && \
        length=(`wc -l {input.samposnohead}`) && \
        if [ ${{length}} -ge 1 ]
        then
        seqtk sample -s2468 {input.R1s} 100000 > {output.R1poshitssub} && \
        seqtk sample -s2468 {input.R2s} 100000 > {output.R2poshitssub} && \
        mmseqs createdb {output.R1poshitssub} {output.R2poshitssub} {params.tempqueryDBfile} --dbtype 2 && \
        mmseqs taxonomy {params.tempqueryDBfile} {params.CO1genome} {params.taxresult} {params.tempqueryDBlocation} \
            -s 1.0 --max-seqs 100 --threads {threads} --search-type 3 --tax-lineage 1 --split-memory-limit 12G --split 2 --remove-tmp-files 1 && \
        mmseqs createtsv {params.tempqueryDBfile} {params.taxresult} {output.taxresulttsv} && \
        mmseqs taxonomyreport {params.CO1genome} {params.taxresult} {output.taxresultkraken} --threads {threads} && \
        rm -r {params.tempqueryDBlocation}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.R1poshitssub} && \
        touch {output.R2poshitssub} && \
        touch {output.taxresulttsv} && \
        touch {output.taxresultkraken}
        """

rule align_all_returned_LSU:
    message:
        """
        microbiome step 2
        Using mmseq2 with LCA based methods to classify taxonomy
        """
    input:
        R1spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_R2s.fastq",
        samposnohead = config["sub_dirs"]["host_detect_LSU"] + "/{sample}_pos_hits_no_head.sam"
    output:
        R1poshitssub = temp(config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_R1_pos_hits_subset.fq"),
        R2poshitssub = temp(config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_R2_pos_hits_subset.fq"),
        taxresulttsv = config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_taxonomy_report_LSU.tsv",
        taxresultkraken = config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_taxonomy_report_LSU_kraken.tsv"
    params:
        LSUgenome = config["LSU_genome_index_mmseq"],
        tempqueryDBfile = config["program_dir"] + "mmeseq/TEMP/LSU/{sample}/{sample}_TEMP_DB",
        tempqueryDBlocation = config["program_dir"] + "mmeseq/TEMP/LSU/{sample}",
        taxresult = config["program_dir"] + "mmeseq/TEMP/LSU/{sample}/{sample}_taxonomy_report"
    log:
        "logs/" + config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}.log"
    threads: 4
    resources:
         mem_mb=24000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.tempqueryDBlocation} ]
        then 
        mkdir -p {params.tempqueryDBlocation} 
        fi && \
        length=$(grep -c -v "^@" {input.samposnohead} 2>/dev/null || echo 0) && \
        if [ ${{length}} -ge 1 ]
        then
        seqtk sample -s2468 {input.R1spos} 100000 > {output.R1poshitssub} && \
        seqtk sample -s2468 {input.R2spos} 100000 > {output.R2poshitssub} && \
        mmseqs createdb {output.R1poshitssub} {output.R2poshitssub} {params.tempqueryDBfile} --dbtype 2 && \
        mmseqs taxonomy {params.tempqueryDBfile} {params.LSUgenome} {params.taxresult} {params.tempqueryDBlocation} \
            -s 1.0 --max-seqs 75 --threads {threads} --search-type 3 --tax-lineage 1 --split-memory-limit 20G --split 2 --remove-tmp-files 1 && \
        mmseqs createtsv {params.tempqueryDBfile} {params.taxresult} {output.taxresulttsv} && \
        mmseqs taxonomyreport {params.LSUgenome} {params.taxresult} {output.taxresultkraken} --threads {threads} && \
        rm -r {params.tempqueryDBlocation}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.R1poshitssub} && \
        touch {output.R2poshitssub} && \
        touch {output.taxresulttsv} && \
        touch {output.taxresultkraken}
        """

rule align_all_returned_SSU:
    message:
        """
        microbiome step 2
        Using mmseq2 with LCA based methods to classify taxonomy
        """
    input:
        R1spos = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R1s.fastq",
        R2spos = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_R2s.fastq",
        samposnohead = config["sub_dirs"]["host_detect_SSU"] + "/{sample}_pos_hits_no_head.sam"
    output:
        R1poshitssub = temp(config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_R1_pos_hits_subset.fq"),
        R2poshitssub = temp(config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_R2_pos_hits_subset.fq"),
        taxresulttsv = config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_taxonomy_report_SSU.tsv",
        taxresultkraken = config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_taxonomy_report_SSU_kraken.tsv"
    params:
        SSUgenome = config["SSU_genome_index_mmseq"],
        tempqueryDBfile = config["program_dir"] + "mmeseq/TEMP/SSU/{sample}/{sample}_TEMP_DB",
        tempqueryDBlocation = config["program_dir"] + "mmeseq/TEMP/SSU/{sample}",
        taxresult = config["program_dir"] + "mmeseq/TEMP/SSU/{sample}/{sample}_taxonomy_report"
    log:
        "logs/" + config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}.log"
    threads: 4
    resources:
         mem_mb=24000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.tempqueryDBlocation} ] 
        then 
        mkdir -p {params.tempqueryDBlocation} 
        fi && \
        length=$(grep -c -v "^@" {input.samposnohead} 2>/dev/null || echo 0) && \
        if [ ${{length}} -ge 1 ] 
        then 
        seqtk sample -s2468 {input.R1spos} 100000 > {output.R1poshitssub} && \
        seqtk sample -s2468 {input.R2spos} 100000 > {output.R2poshitssub} && \
        mmseqs createdb {output.R1poshitssub} {output.R1poshitssub} {params.tempqueryDBfile} --dbtype 2 && \
        mmseqs taxonomy {params.tempqueryDBfile} {params.SSUgenome} {params.taxresult} {params.tempqueryDBlocation} \
            -s 1.0 --max-seqs 65 --threads {threads} --search-type 3 --tax-lineage 1 --split 2 --split-memory-limit 20G --remove-tmp-files 1 && \
        mmseqs createtsv {params.tempqueryDBfile} {params.taxresult} {output.taxresulttsv} && \
        mmseqs taxonomyreport {params.SSUgenome} {params.taxresult} {output.taxresultkraken} --threads {threads} && \
        rm -r {params.tempqueryDBlocation} 
        fi && \
        if [ ${{length}} -eq 0 ] 
        then 
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.R1poshitssub} && \
        touch {output.R2poshitssub} && \
        touch {output.taxresulttsv} && \
        touch {output.taxresultkraken}
        """

rule Generate_figures_CO1_mapping:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from CO1
        """
    input:
        taxresultkraken = config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_taxonomy_report_CO1_kraken.tsv"
    output:
        sampfinished = config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}_top_host_species_additional_stats_CO1.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_LCA_plots_CO1"] + "/",
        nHost = config["Hostdetect"],
        marker = "CO1"
    log:
        "logs/" + config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.taxresultkraken}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/Generate_LCAs_barcode_markers.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        Rscript {config[program_dir]}scripts/emptyLCAs.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        touch {output.sampfinished}
        """

rule Generate_figures_LSU_mapping:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from LSU
        """
    input:
        taxresultkraken = config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_taxonomy_report_LSU_kraken.tsv"
    output:
        sampfinished = config["sub_dirs"]["host_LCA_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_LCA_plots_LSU"] + "/",
        nHost = config["Hostdetect"],
        marker = "LSU"
    log:
        "logs/" + config["sub_dirs"]["host_LCA_plots_LSU"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_LCA_plots_LSU"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.taxresultkraken}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/Generate_LCAs_barcode_markers.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        Rscript {config[program_dir]}scripts/emptyLCAs.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        touch {output.sampfinished}
        """

rule Generate_figures_SSU_mapping:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from SSU
        """
    input:
        taxresultkraken = config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_taxonomy_report_SSU_kraken.tsv"
    output:
        sampfinished = config["sub_dirs"]["host_LCA_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt"
    params:
        samplename = "{sample}",
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["host_LCA_plots_SSU"] + "/",
        nHost = config["Hostdetect"],
        marker = "SSU"
    log:
        "logs/" + config["sub_dirs"]["host_LCA_plots_SSU"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["host_LCA_plots_SSU"] + "/{sample}.txt"
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.taxresultkraken}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/Generate_LCAs_barcode_markers.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        Rscript {config[program_dir]}scripts/emptyLCAs.R \
            --inputtax {input.taxresultkraken} --name {params.samplename} --outputpath {params.wrkdir} \
            --programdir {params.basedir} --marker {params.marker} \
            --host {params.nHost} 
        fi && \
        touch {output.sampfinished}
        """

rule collate_barcodes_results:
    message:
        """
        microbiome step 3
        Generate plots for results and identify top host/bloodmeal species from SSU
        """
    input:
        sampfinishedSSU = config["sub_dirs"]["host_LCA_plots_SSU"] + "/{sample}_top_host_species_additional_stats_SSU.txt",
        sampfinishedLSU = config["sub_dirs"]["host_LCA_plots_LSU"] + "/{sample}_top_host_species_additional_stats_LSU.txt",
        sampfinishedCO1 = config["sub_dirs"]["host_LCA_plots_CO1"] + "/{sample}_top_host_species_additional_stats_CO1.txt"
    output:
        analysisfinished = config["sub_dirs"]["progress_barcodes"] + "/{sample}_completed_fullbarcoding_analysis.txt"
    threads: 1
    shell:
        """
        touch {output.analysisfinished}
        """

####

rule Blast_CO1_contigs:
    message:
        """
        microbiome step 4
        Classifying contig taxonomy of CO1 regions using blastn
        """
    input:
        contigsfile = config["sub_dirs"]["host_contigs_CO1"] + "/{sample}/final.contigs.fa"
    output:
        matchesfile = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["CO1_contig_blast_dir"]
    log:
        "logs/" + config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}.log"
    threads: 1
    resources:
        mem_mb=6000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        length=(`wc -l {input.contigsfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        blastn -query {input.contigsfile} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 1 \
            -word_size 22 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -out {output.matchesfile} \
            2> {log} && \
        touch {output.matchesfile}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.matchesfile}
        """

rule Blast_LSU_contigs:
    message:
        """
        microbiome step 4
        Classifying contig taxonomy of LSU regions using blastn 
        """
    input:
        contigsfile = config["sub_dirs"]["host_contigs_LSU"] + "/{sample}/final.contigs.fa"
    output:
        matchesfile = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["LSU_contig_blast_dir"]
    log:
        "logs/" + config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}.log"
    threads: 1
    resources:
        mem_mb=4000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        length=(`wc -l {input.contigsfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        blastn -query {input.contigsfile} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 1 \
            -max_hsps 1 \
            -word_size 22 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -out {output.matchesfile} \
            2> {log} && \
        touch {output.matchesfile}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.matchesfile}
        """

rule Blast_SSU_contigs:
    message:
        """
        microbiome step 4
        Classifying contig taxonomy of SSU regions using blastn 
        """
    input:
        contigsfile = config["sub_dirs"]["host_contigs_SSU"] + "/{sample}/final.contigs.fa"
    output:
        matchesfile = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["SSU_contig_blast_dir"]
    log:
        "logs/" + config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}.log"
    threads: 1
    resources:
        mem_mb=4000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        length=(`wc -l {input.contigsfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        blastn -query {input.contigsfile} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 1 \
            -word_size 22 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -out {output.matchesfile} \
            2> {log} && \
        touch {output.matchesfile}
        fi && \
        if [ ${{length}} -eq 0 ]
        then
        echo " No hits were found for this barcode marker region. This may indicate there is insufficient information to determine host from just barcode markers and manual inspection of results is recommended. Pipeline will still attempt host identification but consider the likely accuracy of the chosen host "
        fi && \
        touch {output.matchesfile}
        """

rule analyse_blastn_CO1:
    message:
        """
        DNA classification step 2
        generate taxonomies from each blastn hit for {wildcards.sample} using R and taxonomizr package
        """
    input:
        matchesfile = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    output:
        sampfinished = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_finished.txt",
        unassignedDNAcontigs = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_matches_unassigned.fa",
        contigsassignedfile = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_Contigsallinformationassignment.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["CO1_contig_blast_dir"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["CO1_contig_blast_dir"] + "/",
        contigname = "{sample}_matches_unassigned.fa"
    log:
        "logs/" + config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_assignments.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_assignments.txt"
    resources:
        mem_mb=2000
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.matchesfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_best_hits_microbiome_HPC.R \
            --inputblastn {input.matchesfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.unassignedDNAcontigs} --programdir {params.basedir} --savdir {params.wrkdir} \
            --abundances {params.abundances} --savcontig {params.contigname}
        fi && \
        touch {output.sampfinished} && \
        touch {output.contigsassignedfile} && \
        touch {output.unassignedDNAcontigs}
        """

rule analyse_blastn_LSU:
    message:
        """
        DNA classification step 2
        generate taxonomies from each blastn hit for {wildcards.sample} using R and taxonomizr package
        """
    input:
        matchesfile = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    output:
        sampfinished = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_finished.txt",
        unassignedDNAcontigs = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_matches_unassigned.fa",
        contigsassignedfile = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_Contigsallinformationassignment.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["LSU_contig_blast_dir"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["LSU_contig_blast_dir"],
        contigname = "{sample}_matches_unassigned.fa"
    log:
        "logs/" + config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_assignments.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_assignments.txt"
    resources:
        mem_mb=2000
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.matchesfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_best_hits_microbiome_HPC.R \
            --inputblastn {input.matchesfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.unassignedDNAcontigs} --programdir {params.basedir} --savdir {params.wrkdir} \
            --abundances {params.abundances} --savcontig {params.contigname}
        fi && \
        touch {output.sampfinished} && \
        touch {output.contigsassignedfile} && \
        touch {output.unassignedDNAcontigs}
        """

rule analyse_blastn_SSU:
    message:
        """
        DNA classification step 2
        generate taxonomies from each blastn hit for {wildcards.sample} using R and taxonomizr package
        """
    input:
        matchesfile = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_matches_nucleotide.m8"
    output:
        sampfinished = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_finished.txt",
        unassignedDNAcontigs = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_matches_unassigned.fa",
        contigsassignedfile = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_Contigsallinformationassignment.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["SSU_contig_blast_dir"],
        basedir = config["program_dir"],
        wrkdir = config["sub_dirs"]["SSU_contig_blast_dir"],
        contigname = "{sample}_matches_unassigned.fa"
    log:
        "logs/" + config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_assignments.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_assignments.txt"
    resources:
        mem_mb=2000
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.matchesfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_best_hits_microbiome_HPC.R \
            --inputblastn {input.matchesfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --output {output.unassignedDNAcontigs} --programdir {params.basedir} --savdir {params.wrkdir} \
            --abundances {params.abundances} --savcontig {params.contigname}
        fi && \
        touch {output.sampfinished} && \
        touch {output.contigsassignedfile} && \
        touch {output.unassignedDNAcontigs}
        """


rule collate_barcodes_results_contigs:
    message:
        """
        microbiome step 3
        Collating results for microbiome analysis contigs generated
        """
    input:
        sampfinishedC01 = config["sub_dirs"]["CO1_contig_blast_dir"] + "/{sample}_finished.txt",
        sampfinishedLSU = config["sub_dirs"]["LSU_contig_blast_dir"] + "/{sample}_finished.txt",
        sampfinishedSSU = config["sub_dirs"]["SSU_contig_blast_dir"] + "/{sample}_finished.txt"
    output:
        analysisfinished = config["sub_dirs"]["progress_barcodes"] + "/{sample}_completed_fullbarcoding_analysis_contigs.txt"
    threads: 1
    shell:
        """
        touch {output.analysisfinished}
        """