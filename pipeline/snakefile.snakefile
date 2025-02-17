import os, sys

configfile: "config.yaml"
"""
# Call the relevant variables which control which rules will be used in the analysis
"""

hostdetection = config["inspecthost"]

DNA_analysis = config["DNA_assign_blastn"]
Blastmethod = config["Blast_method"]

Blastn_false_positive = config["Blastn_viral_contig_false_positive_check"]

Advanced_microbiome = config["Microbiome_classification"]

Advanced_viral_genome_build = config["Viral_genome_build"]
Viral_trees = config["Viral_genome_tree_building"]
genomad_run = config["Genomad_detect"]

analyse_raws = config["run_raws"]
diamondraw = config["dodiamond_blast_raws"]
divergent_reads_and_contigs = config["Divergent_reads_and_contigs_search"]

summary = config["generate_results_summary"]

"""
# Define the base rule all plus the rule all's for each specific branch of the analysis (e.g., host mapping).
# This is the proper final one
"""
rule_all_input_list=expand(config["sub_dirs"]["raws_to_contigs"] + "/{sample}_hits.sam", sample=config["samples"])
rule_all_advanced_microbiome=expand(config["sub_dirs"]["progress_barcodes"] + "/{sample}_completed_fullbarcoding_analysis.txt", sample=config["samples"])
rule_all_DNA_analysis=expand(config["sub_dirs"]["contigs_assigned_nucl_abundances"] + "/{sample}_finished.txt", sample=config["samples"])
rule_all_host_removal=expand(config["sub_dirs"]["raws_to_host_contigs"] + "/{sample}_host_contig_hits.sam", sample=config["samples"])
rule_all_summary=expand(config["sub_dirs"]["finished"] + "/{sample}_finished",sample=config["samples"])
rule_all_advanced_microbiome_contigs=expand(config["sub_dirs"]["progress_barcodes"] + "/{sample}_completed_fullbarcoding_analysis_contigs.txt", sample=config["samples"])
rule_all_genomad=expand(config["sub_dirs"]["Genomad_metabat_sorted"] + "/{sample}/finished.txt", sample=config["samples"])

rule_all_diamond_raws=expand(config["sub_dirs"]["compiled_summary"] + "/{sample}/finished_extracting_reads.txt",sample=config["samples"])
rule_all_diamond_raws_contigs_diverged=expand(config["sub_dirs"]["compiled_summary"] + "/{sample}/finished_extracting_reads_diverged.txt",sample=config["samples"])

rule_all_viral_genome_building=expand(config["sub_dirs"]["Viral_genomes_present_refined_genomes"] + "/{sample}/finished2.txt", sample=config["samples"])
rule_all_viral_genome_trees=expand(config["sub_dirs"]["Viral_genomes_present_iqtree"] + "/{sample}/finished.txt", sample=config["samples"])
rule_all_blastn_false_positive_check=expand(config["sub_dirs"]["finished"] + "/{sample}_finished_blastn_false_positive_check", sample=config["samples"])

"""
# define the base path to all rules plus the main rules always present (those in the base rule all)
"""
rules_dir = os.path.join(os.path.expanduser(config["program_dir"]), "rules")
include: os.path.join(rules_dir, "preprocessing.smk")
include: os.path.join(rules_dir, "assembly_analysis.smk")
"""
# append/extend depending on the size of the lists to be added
"""
if config["Host_filter"] == 'yes':
    rule_all_input_list.extend(rule_all_host_removal)
    include: os.path.join(rules_dir, "identify_host_sp.smk")
    include: os.path.join(rules_dir, "Microbiome_advanced.smk")
    print("The host species will be identified and filtered out before main analysis")
# remove mmseq2 from DNA analysis stage
if DNA_analysis == 'yes':

    rule_all_input_list.extend(rule_all_DNA_analysis)
    
    if Blastmethod == 'mmseq2':
        include: os.path.join(rules_dir, "mmseq2_DNA_classify_blastn.smk")
        print("nucleotide blast will be performed on contigs using mmseq2")

    elif Blastmethod == 'blast':
        include: os.path.join(rules_dir, "DNA_classify.smk")
        print("nucleotide blast will be performed on contigs using blastn")

else:
    print("unassigned contigs will not be nucleotide blasted")

if Advanced_microbiome == 'yes':

    rule_all_input_list.extend(rule_all_advanced_microbiome)
    rule_all_input_list.extend(rule_all_advanced_microbiome_contigs)
    include: os.path.join(rules_dir, "Microbiome_advanced.smk")
    print("advanced rRNA microbiome classification will be performed")

else:
    print("advanced LCA analysis of the microbiome LSU and SSU regions will not be undertaken")


if Advanced_viral_genome_build == "yes":

    rule_all_input_list.extend(rule_all_viral_genome_building)
    include: os.path.join(rules_dir, "Generate_virus_genomes.smk")
    print("Consensus genomes of detected viruses will be generated")

else:
    print("estimation of genome coverage and of Viral genomes will not be undertaken")

if Viral_trees == "yes":

    rule_all_input_list.extend(rule_all_viral_genome_trees)
    include: os.path.join(rules_dir, "Generate_trees_virus_genomes.smk")
    print("Consensus genomes of detected viruses will aligned to make trees")

else:
    print("Generating phylogenetic trees of viruses not undertaken")

if Blastn_false_positive == 'yes':

    rule_all_input_list.extend(rule_all_blastn_false_positive_check)
    include: os.path.join(rules_dir, "Blastn_false_positive_check.smk")
    print("Viral contigs identified through Diamond BLASTx will be checked against BLASTn for better matches to minimise false positives")

else:
    print("Secondary viral confirmation through BLASTn will not be undertaken")


if genomad_run == 'yes':

    rule_all_input_list.extend(rule_all_genomad)
    include: os.path.join(rules_dir, "Genomad_viral_assignments.smk")
    print("Genomad will be run and viral gene regions likely from the same species will be collated")

else:
    print("Neural network identification of likely viral fragments by genomad will not be undertaken")

if summary == 'yes':

    rule_all_input_list.extend(rule_all_summary)
    include: os.path.join(rules_dir, "summarise_results.smk")
    print("Summary results will be generated")

else:
    print("no summary results will be generated")

if analyse_raws == 'yes':

    include: os.path.join(rules_dir, "analyse_raws.smk")
    print("the raw sequences not assigned to a generated contig will be analysed using the set programs")

else:
    print("reads not assigned to contigs will not be analysed further")

if diamondraw  == 'yes':
    rule_all_input_list.extend(rule_all_diamond_raws)
    print("Raw sequences will be analysed using Diamond protein blast")


if divergent_reads_and_contigs  == 'yes':
    rule_all_input_list.extend(rule_all_diamond_raws_contigs_diverged)
    include: os.path.join(rules_dir, "Diverged_read_detection_viruses.smk")
    print("Taxids of detected viruses will be used to perform an more sensitive viral search on reads and contigs for previously detected viruses")


rule all:
    input:
        data = rule_all_input_list
    output:
        finished = config["program_dir"] + "All_samples_finished_analysis.txt"
    shell:
        """
        touch {output.finished}
        """
