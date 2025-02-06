
def get_raw_fastqs(wildcards):
    return config["samples"][wildcards.sample]


def Assembly_used_initfilter(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_host_removed.fa"
        ])
    if config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["contig_dir_either"] + "/{sample}_contigs.fa"
        ])


rule Prep_viral_contigs:
    message:
        """
        Subsetting contigs to just viral assigned contigs for secondary confirmation with blastn 
        subsetting contigs for {wildcards.sample}
        """
    input:
        contigsassignedtab = config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_Contigsallinformationassignment.txt",
        readslist2 = config["sub_dirs"]["Kraken_viral_ids_contigs"] + "/{sample}_contigs_to_virus.txt",
        contigfile = Assembly_used
    output:
        contigslist = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_viral_assigned_temp1.tsv",
        contigslist2 = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_viral_assigned_contig_names.tsv",
        virus_contigs = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_viral_assigned_contigs.fasta"
    threads: 4
    resources:
        mem_mb=4000
    shell:
        """
        awk -F '\t' '$12 == "Viruses" {{print $1}}' "{input.contigsassignedtab}" > {output.contigslist} && \
        cat {output.contigslist} {input.readslist2} >> {output.contigslist2}
        if [ -n "{output.contigslist}" ]; then
        seqkit grep -f {output.contigslist2} {input.contigfile} > {output.virus_contigs}
        fi && \
        touch {output.virus_contigs} && \
        touch {output.contigslist}
        """

rule false_positive_check_blastn:
    message:
        """
        DNA classification step 1 (starts after analysis step 3) 
        Blast unassigned contigs {wildcards.sample} blastn
        """
    input:
        virus_contigs = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_viral_assigned_contigs.fasta"
    output:
        blastfile = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_matches_nucleotide.m8"
    params:
        blastdb= config["blast_nucleotide_database"],
        blastnfolder = config["sub_dirs"]["contigs_nucl_false_pos_check"]
    log:
        "logs/" + config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}.log"
    threads: 4
    resources:
        mem_mb=32000
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}.txt"
    shell:
        """
        if [ ! -d {params.blastnfolder} ]
        then
            mkdir {params.blastnfolder}
        fi && \
        lengthunassigned=(`wc -l {input.virus_contigs}`) && \
        if [ ${{lengthunassigned}} -ge 2 ]
        then
        blastn -query {input.virus_contigs} \
            -db {params.blastdb} \
            -evalue 0.00001 \
            -max_target_seqs 4 \
            -max_hsps 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp' \
            -num_threads {threads} \
            -word_size 24 \
            -out {output.blastfile} \
            2> {log}
        fi && \
        touch {output.blastfile}
        """

rule analyse_blastn_hits_false_positive:
    message:
        """
        DNA classification step 2
        generate taxonomies from each Diamond hit for {wildcards.sample} raw reads using R and taxonomizr package
        """
    input:
        blastfile = config["sub_dirs"]["contigs_nucl_false_pos_check"] + "/{sample}_matches_nucleotide.m8"
    output:
        sampfinished = config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/{sample}_finished.txt",
        contigsassignedfile = config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/abundances/{sample}_Contigsallinformationassignment.txt"
    params:
        samplename = "{sample}",
        namenodedatabase = config["Accession_allnamenode"],
        abundances = config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/abundances/",
        basedir = config["program_dir"],
    log:
        "logs/" + config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/{sample}.txt"
    resources:
        mem_mb=2000
    conda: "Rdataplotting"
    threads: 1
    shell:
        """
        length=(`wc -l {input.blastfile}`) && \
        if [ ${{length}} -ge 1 ]
        then
        Rscript {config[program_dir]}scripts/extract_blastn_false_positive_virus.R \
            --inputblastn {input.blastfile} --name {params.samplename} \
            --threads {threads} --Accnode {params.namenodedatabase} --programdir {params.basedir} \
            --abundances {params.abundances}
        fi && \
        touch {output.sampfinished} && \
        touch {output.contigsassignedfile}
        """

rule Add_blastn_false_positive_check:
    message:
        """
        Optional summary step
        Adding in blastn false positive check results
        """
    input:
        summaryRdatafile = config["sub_dirs"]["Summary_results"] + "/{sample}_gather_summary_files_R_environment.Rdata",
        contigsassignedfile = config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/abundances/{sample}_Contigsallinformationassignment.txt",
    output:
        Revised_table_viruses = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly_including_blastn_false_positive_check.txt"
    params:
        inputpath = config["sub_dirs"]["Summary_results"] + "/",
        outputpth = config["sub_dirs"]["Summary_results"] + "/",
        basedir = config["program_dir"],
        sampname = "{sample}"
    conda: "Rdataplotting"
    threads: 1
    resources:
         mem_mb=2000
    shell:
        """
        Rscript {config[program_dir]}scripts/generate_summary_tables_gather_results_false_positive_check.R \
            --inputRenv {input.summaryRdatafile} \
            --inputblastn_results {input.contigsassignedfile} \
            --outputpath {params.outputpth} --programdir {params.basedir} \
            --inputpath {params.inputpath} --samplename {params.sampname} && \
        touch {output.Revised_table_viruses}
        """

rule Compile_shared_graphs_added_blastn_false_positive:
    message:
        """
        Optional summary step 2
        Adding in false positive blastn check to final combined plots and tables across samples
        """
    input:
        summaryallreadsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}Summary_assignment_reads_for_plot_generation.txt",sample=config["samples"]),
        summaryreadsfilterfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summary_reads_filtering.txt",sample=config["samples"]),
        summarycontigsstatsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontigs_assembly_values.txt",sample=config["samples"]),
        summaryrawdiamondhitsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt",sample=config["samples"]),
        Revised_table_viruses_files = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly_including_blastn_false_positive_check.txt",sample=config["samples"])
    output:
        completeRenv = config["sub_dirs"]["Summary_results2"] + "/gather_summary_files_R_environment_false_positive_blastncheck.Rdata"
    params:
        inputpath = config["sub_dirs"]["Summary_results"] + "/",
        dohostdetect = config["Host_filter"],
        dodiamondraws = config["dodiamond_blast_raws"],
        outputpth = config["sub_dirs"]["Summary_results2"] + "/",
        basedir = config["program_dir"]
    log:
        "logs/" + config["sub_dirs"]["Summary_results2"] + "/summarygraphs.log"
    conda: "Rdataplotting"
    threads: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Summary_results2"] + "/summarygraphs.txt"
    resources:
         mem_mb=2000
    shell:
        """
        Rscript {config[program_dir]}scripts/Generate_combined_Plots_and_tables_blastn_false_positive_check.R \
            --dohostdetect {params.dohostdetect} --dodiamondraws {params.dodiamondraws} \
            --outputpath {params.outputpth} --programdir {params.basedir} \
            --inputpath {params.inputpath} && \
        touch {output.completeRenv}
        """



rule finished_blastn_false_positive:
    input:
        completeRenv = config["sub_dirs"]["Summary_results2"] + "/gather_summary_files_R_environment_false_positive_blastncheck.Rdata",
        summaryrawdiamondhitsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt"
    output:
        finished = config["sub_dirs"]["finished"] + "/{sample}_finished_blastn_false_positive_check"
    shell:
        """
        touch {output.finished}
        """


