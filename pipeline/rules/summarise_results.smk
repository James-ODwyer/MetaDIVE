
def Assembly_used_log(wildcards):
    if config["Assembly_choice"] == "Trinity":
        return([
            config["sub_dirs"]["contig_dir_trinity"] + "/{sample}_trinity.Trinity.fasta"
        ])
    elif config["Assembly_choice"] == "Megahit":
        return([
            "logs/" + config["sub_dirs"]["contig_dir_megahit"] + "/{sample}.log"
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

# Fix below still. should be for the reads aligned to host contigs not the host contigs to host genome
def hostremovedreads(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            "logs/" + config["sub_dirs"]["host_species_genomes"] + "/aligning/{sample}.txt"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])

def hostcontigsrawsalign(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            "logs/" + config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}.log"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])


def hostremovedcontigsfa(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}_contigs_aligned_host.fa"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])

def hostremovedcontigssam(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            "logs/" + config["sub_dirs"]["contig_dir_host_rem"] + "/{sample}.txt"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])

#config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_host_contig_hits.sam"

def hostspecies(wildcards):
    if config["Host_filter"] == 'yes':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/{sample}_top_host_species_overall.txt"
        ])
    elif config["Host_filter"] == 'no':
        return([
            config["sub_dirs"]["host_species_genomes"] + "/No_host_identification_undertaken_{sample}.txt"
        ])

def diamondrawskingdoms(wildcards):
    if config["dodiamond_blast_raws"] == 'yes':
        return([
            config["sub_dirs"]["contigs_assigned_raw"] + "/{sample}_contigs_file.txt"
        ])
    elif config["dodiamond_blast_raws"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def diamondrawsalignments(wildcards):
    if config["dodiamond_blast_raws"] == 'yes':
        return([
            "logs/" + config["sub_dirs"]["raws_to_contigs"] + "/{sample}.log"
        ])
    elif config["dodiamond_blast_raws"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def krakenraws(wildcards):
    if config["doviral_raws_classification"] == 'yes':
        return([
            config["sub_dirs"]["Kraken_Vir_raws"] + "/{sample}_kraken_report_minimiserdata.tsv"
        ])
    elif config["doviral_raws_classification"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def DNA_assign_blast_contigs(wildcards):
    if config["DNA_assign_blastn"] == 'yes':
        return([
            config["sub_dirs"]["contigs_assigned_nucl"] + "/abundances/{sample}_Contigsallinformationassignment.txt"
        ])
    elif config["DNA_assign_blastn"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def MicrobiomeCO1(wildcards):
    if config["Microbiome_classification"] == 'yes':
        return([
            config["sub_dirs"]["host_assignLCA_CO1"] + "/{sample}_taxonomy_report_CO1_kraken.tsv"
        ])
    elif config["Microbiome_classification"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def MicrobiomeLSU(wildcards):
    if config["Microbiome_classification"] == 'yes':
        return([
            config["sub_dirs"]["host_assignLCA_LSU"] + "/{sample}_taxonomy_report_LSU_kraken.tsv"
        ])
    elif config["Microbiome_classification"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def MicrobiomeSSU(wildcards):
    if config["Microbiome_classification"] == 'yes':
        return([
            config["sub_dirs"]["host_assignLCA_SSU"] + "/{sample}_taxonomy_report_SSU_kraken.tsv"
        ])
    elif config["Microbiome_classification"] == 'no':
        return([
            "PROGRESS_STATUS/No_analysis_undertaken_{sample}.txt"
        ])

def Diamondresultscontigsviruscheck(wildcards):
    if config["Blastn_viral_contig_false_positive_check"] == 'yes':
        return([
            config["sub_dirs"]["contigs_nucl_false_pos_check_counts"] + "/abundances/{sample}_Contigsallinformationassignment.txt"
        ])
    elif config["Blastn_viral_contig_false_positive_check"] == 'no':
        return([
            config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_Contigsallinformationassignment.txt"
        ])

rule summarise_all_results:
    message:
        """
        Summary step
        Collecting results from all previous analyses
        """
    input:
        Fastplog = "logs/" + config["sub_dirs"]["trim_dir"] + "/{sample}.log",
        Phixlog = "logs/" + config["sub_dirs"]["PhiX_dir"] + "/{sample}.log",
        CO1log = "logs/" + config["sub_dirs"]["CO1_dir"] + "/{sample}.log",
        LSUlog = "logs/" + config["sub_dirs"]["LSU_dir"] + "/{sample}.log",
        SSUlog = "logs/" + config["sub_dirs"]["SSU_dir"] + "/{sample}.log",
        CO1assignments = MicrobiomeCO1,
        LSUassignments = MicrobiomeLSU,
        SSUassignments = MicrobiomeSSU,
        Assembly_log = Assembly_used_log,
        hostsp = hostdetected,
        Diamondresultscontigsviruscheckfalsepos = Diamondresultscontigsviruscheck,
        diamondresultscontigs = config["sub_dirs"]["contigs_assigned"] + "/abundances/{sample}_Contigsallinformationassignment.txt",
        DNA_contigs_results = DNA_assign_blast_contigs,
        hostalignmentrates = hostremovedreads,
        host_removed_contigsfa = hostremovedcontigsfa,
        hostcontigsrawsbowtie = hostcontigsrawsalign,
        readsassignedhostcontigs = hostremovedcontigssam,
        diamondrawskingdoms = diamondrawskingdoms,
        diamondrawsalign = diamondrawsalignments,
        krakenrawstats = krakenraws,
        raw_log = "logs/" + config["sub_dirs"]["raws_to_contigs"] + "/{sample}.log",
        hostspecieschosen = hostspecies,
        samrawscontigs = config["sub_dirs"]["raws_to_contigs"] + "/{sample}_hits_3col.sam"
    output:
        summaryfile = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assembly.txt",
        summaryallreadsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}Summary_assignment_reads_for_plot_generation.txt",
        summaryreadsfilterfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summary_reads_filtering.txt",
        summarycontigsstatsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontigs_assembly_values.txt",
        summaryrawdiamondhitsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt",
        summaryrawdiamondhitsreadnames = config["sub_dirs"]["Summary_results"] + "/{sample}_raw_read_names_to_virus.txt",
        summaryRdatafile = config["sub_dirs"]["Summary_results"] + "/{sample}_gather_summary_files_R_environment.Rdata"
    params:
        sampname = "{sample}",
        dohostdetect = config["Host_filter"],
        dodiamondraws = config["dodiamond_blast_raws"],
        dokrakenraws = config["doviral_raws_classification"],
        doDNAblastn = config["DNA_assign_blastn"],
        domicrobiome = config["Microbiome_classification"],
        outputpth = config["sub_dirs"]["Summary_results"] + "/",
        basedir = config["program_dir"],
        namenodedatabase = config["Accession_allnamenode"],
        Assemblychoice = config["Assembly_choice"],
        contigfalseposchoice = config["Blastn_viral_contig_false_positive_check"],
        doonlyblastncontigs = config["Final_contigs_returned"],
        nohostgenomemarker = "PROGRESS_STATUS/No_host_genome_was_found_{sample}.txt"
    log:
        "logs/" + config["sub_dirs"]["Summary_results"] + "/{sample}.log"
    conda: "Rdataplotting"
    threads: 1
    benchmark:
        "benchmarks/" + config["sub_dirs"]["Summary_results"] + "/{sample}.txt"
    resources:
         mem_mb=24000
    shell:
        """
        if [ ! -f {params.nohostgenomemarker} ]
        then
        Rscript {config[program_dir]}scripts/gather_results.R \
            --fastplog {input.Fastplog} --phiXlog {input.Phixlog} --CO1_bowtie {input.CO1log} --LSU_bowtie {input.LSUlog} --SSU_bowtie {input.SSUlog} \
            --CO1microbiome {input.CO1assignments} --LSUmicrobiome {input.LSUassignments} --SSUmicrobiome {input.SSUassignments} \
            --Assembly_log {input.Assembly_log} --hostsp {input.hostsp} --hostaligngenome {input.hostalignmentrates} --Assemblyused {params.Assemblychoice} \
            --diamondrawskingdoms {input.diamondrawskingdoms} --diamondrawslog {input.diamondrawsalign} --krakenraws {input.krakenrawstats} --rawtocontigsassignlog {input.raw_log} \
            --samplename {params.sampname} --dohostdetect {params.dohostdetect} --dodiamondraws {params.dodiamondraws} --identifiedhost {input.hostspecieschosen} \
            --dokrakenraws {params.dokrakenraws} --doblastn {params.doDNAblastn} --domicrobiome {params.domicrobiome} --dohostcontigsrawsalign {input.hostcontigsrawsbowtie} \
            --outputpath {params.outputpth} --Log {log} --programdir {params.basedir} --rawreadssamcontigs {input.samrawscontigs} --hostalignmentreads {input.hostalignmentrates} \
            --Diamondtab {input.diamondresultscontigs} --Blastntab {input.DNA_contigs_results} --rawssam {input.samrawscontigs} --hostremovedcontigsfa {input.host_removed_contigsfa} --hostremovedcontigssam {input.readsassignedhostcontigs} --hostgenomefailed yes \
            --Accnode {params.namenodedatabase} --dofalseposcontigschoice {params.contigfalseposchoice} --falseposcontigsrem {input.Diamondresultscontigsviruscheckfalsepos} --doblastnassignmentsonly {params.doonlyblastncontigs} && \
        touch {output.summaryfile} && \
        touch {output.summaryallreadsfiles} && \
        touch {output.summaryreadsfilterfiles} && \
        touch {output.summaryrawdiamondhitsfiles} && \
        touch {output.summaryrawdiamondhitsreadnames} && \
        touch {output.summaryRdatafile}
        fi && \
        if [ -f {params.nohostgenomemarker} ]
        then
        Rscript {config[program_dir]}scripts/gather_results.R \
            --fastplog {input.Fastplog} --phiXlog {input.Phixlog} --CO1_bowtie {input.CO1log} --LSU_bowtie {input.LSUlog} --SSU_bowtie {input.SSUlog} \
            --CO1microbiome {input.CO1assignments} --LSUmicrobiome {input.LSUassignments} --SSUmicrobiome {input.SSUassignments} \
            --Assembly_log {input.Assembly_log} --hostsp {input.hostsp} --hostaligngenome {input.hostalignmentrates} --Assemblyused {params.Assemblychoice} \
            --diamondrawskingdoms {input.diamondrawskingdoms} --diamondrawslog {input.diamondrawsalign} --krakenraws {input.krakenrawstats} --rawtocontigsassignlog {input.raw_log} \
            --samplename {params.sampname} --dohostdetect {params.dohostdetect} --dodiamondraws {params.dodiamondraws} --identifiedhost {input.hostspecieschosen} \
            --dokrakenraws {params.dokrakenraws} --doblastn {params.doDNAblastn} --domicrobiome {params.domicrobiome} --dohostcontigsrawsalign {input.hostcontigsrawsbowtie} \
            --outputpath {params.outputpth} --Log {log} --programdir {params.basedir} --rawreadssamcontigs {input.samrawscontigs} --hostalignmentreads {input.hostalignmentrates} \
            --Diamondtab {input.diamondresultscontigs} --Blastntab {input.DNA_contigs_results} --rawssam {input.samrawscontigs} --hostremovedcontigsfa {input.host_removed_contigsfa} --hostremovedcontigssam {input.readsassignedhostcontigs} --hostgenomefailed yes \
            --Accnode {params.namenodedatabase} --dofalseposcontigschoice {params.contigfalseposchoice} --falseposcontigsrem {input.Diamondresultscontigsviruscheckfalsepos} --doblastnassignmentsonly {params.doonlyblastncontigs} && \
        touch {output.summaryfile} && \
        touch {output.summaryallreadsfiles} && \
        touch {output.summaryreadsfilterfiles} && \
        touch {output.summaryrawdiamondhitsfiles} && \
        touch {output.summaryrawdiamondhitsreadnames} && \
        touch {output.summaryRdatafile}
        fi
        """

rule Compile_shared_graphs:
    message:
        """
        Summary step 2
        Generating final combined plots and tables across samples
        """
    input:
        summaryallreadsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}Summary_assignment_reads_for_plot_generation.txt",sample=config["samples"]),
        summaryreadsfilterfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summary_reads_filtering.txt",sample=config["samples"]),
        summarycontigsstatsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontigs_assembly_values.txt",sample=config["samples"]),
        summaryrawdiamondhitsfiles = expand(config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt",sample=config["samples"])
    output:
        completeRenv = config["sub_dirs"]["Summary_results2"] + "/gather_summary_files_R_environment.Rdata"
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
         mem_mb=12000
    shell:
        """
        Rscript {config[program_dir]}scripts/Generate_combined_Plots_and_tables.R \
            --dohostdetect {params.dohostdetect} --dodiamondraws {params.dodiamondraws} \
            --outputpath {params.outputpth} --programdir {params.basedir} \
            --inputpath {params.inputpath} && \
        touch {output.completeRenv}
        """

rule finished:
    input:
        completeRenv = config["sub_dirs"]["Summary_results2"] + "/gather_summary_files_R_environment.Rdata",
        summaryrawdiamondhitsfiles = config["sub_dirs"]["Summary_results"] + "/{sample}_summarycontighits_assigned_assembly.txt"
    output:
        finished = config["sub_dirs"]["finished"] + "/{sample}_finished"
    log:
        "logs/" + config["sub_dirs"]["finished"] + "/{sample}_finished.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["finished"] + "/{sample}_finished.txt"
    shell:
        """
        touch {output.finished}
        """
       
