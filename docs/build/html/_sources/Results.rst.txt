Interpreting results
====================

Results folders
---------------

MetaDIVE will generate slightly different results based on the modules chosen but the main results of all runs will be contained in a results folder titled **Summarised_results**

This folder will contain several subfolders depending on the modules used. Here we will break down the results of each folder and ways to interpret them

.. image:: images/Summarised_results_folder.png
   :alt: MetaDIVE Summary results
   :width: 400px
   :align: center


This main summary results folder is located in the pipeline directory after analysis is complete.


/MetaDIVE/pipeline/Summarised_results



Summary_results_combined_figures_and_tables_contigs
---------------------------------------------------

The following folder contains a breakdown of all summary graphs and tables from MetaDIVE. 


.. image:: images/summary_results_and_tables.png
   :alt: MetaDIVE Summary results
   :width: 400px
   :align: center


All figures in MetaDIVE are saved as both png and pdf files for easy reading.

The key files to be aware of are 


Top_viral_hits_combined_raws_plus_contigs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Top_viral_hits_combined_raws_plus_contigs** (txt or png or pdf) are only present if Single read analysis module is run.
These files contain the most comprehensive summaries of all MetaDIVE analyses for viral assignments. These files 

The txt file of this contains the following column headers 

**reads_assigned_through_contigs**

**superkingdom, phylum class, order,family,genus,species,subspecies** (Taxonomy of closest viral match)

**taxid** (taxID of closest viral match)

**average_percent_ident, min_percent_ident, max_percent_ident** (The average, min and max pairwise identity match of the contig in your sample to the closest matching reference)

**length** (the average length of contigs from your sample which matched to the closest matching reference

**number_contigs_assigned** (The number of contigs which assigned to this viral species)

**mean_complexity_of_contigs** (The average compexity of the DNA sequence of the generated contig which provides insights into how redundant the sequence is)

**contigs_assigned** (binary yes/no for whether contigs were assigned to this species at all, can be ignored)	

**additional_raw_reads_assigned** (binary yes/no for whether single reads were assigned to this species at all, can be ignored)

**number_of_raw_reads_assigned** (how many single reads were assigned to this reference species)
	
**mean_complexity_of_raw_reads** (the average compexity of the DNA sequence of the reads assigned which provides insights into how redundant the sequences are)

**mean_identity_of_raw_reads** (the average pairwise identity match of the reads matching to the the closest matching reference)

**mean_aligned_length_of_raw_reads** (The average alignment length of the reads to this reference viral species from BLASTn)

**total_reads_assigned** (the total number of reads assigned through contigs [all which assigned to that viral species] and reads [only those which assigned through BLASTx and BLASTn] )

**classification_level_estimate** (informal estimate of the likelihood of the viral species assigned being the same species as the reference, this is still a work in progress)

**confidence_of_assignment_as_virus** (informal estimate of the likelihood of the viral species identified being in the data at all vs a false positive, this is still a work in progress)

**Sample** (Your sample name)


From this results file you can see a comprehensive summary of what was detected and how likely what was assigned as the closest match is to be true. 


Summarybarplot_contigs(bacterial,Eukaryote) and Combined_samples_top_hits_to(bacterial,Eukaryote)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files and figures contain the top 10 identified bacterial and Eukaryote species (excluding host if host identification module was used). 
These are similar but with fewer columns than the above example of Top_viral_hits_combined_raws_plus_contigs. Caution should be used when reading these tables as bacterial and Eukaryote assignments are prone to false positive assignments
through blastx. 


Summarybarplot_(bacterial,Eukaryote,virus)_families
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files and figures are the like the above summarybarplot contigs results except with just single read results from Diamond BLASTx. This is done before Single Read analysis and so no quality 
checks are completed yet. These files are only useful if you do not want to run the full Single Read analysis module but want a rough idea of what the rest of the reads are composed of.

Combined_samples_top_hits_to_viral_species
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are multiple similar tables here which serve slightly different functions. 

**Combined_samples_top_hits_to_Viral_species.txt** (pairwise count matrix of reads assigned to top 10 viruses across all samples from only contigs)

**Combined_samples_top_hits_to_Viral_species_long_table_format.txt** (As above but instead of a matrix is in a longer format with additional information of pairwise identity of alignments and alignment lengths)

**Combined_samples_top_hits_to_Viral_species_blastn_false_positive_check.txt** and **Combined_samples_top_hits_to_Viral_species_long_table_format_blastn_false_positive_check.txt** (like the above two tables, but with added information on comparsions between BLASTx and BLASTn results to 
help tease apart unusual results.


Combined_summary_reads_filtering.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary of all the reads assignments as part of the core module e.g, how many reads passed QC. were assigned as CO1, to host, to SSU, to bacteria, Viruses etc.

paired text file to **summary_all_reads_assigned_and_filtered.(png/pdf)** which shows this information in figure format

Similar also to **Combined_summary_contigs_assiging.txt** which has the same information but focussed on contigs, contig sizes, contig numbers, numbers assigned through 
BLASTx, BLASTn or unassigned.


gather_summary_files_R_environment and gather_summary_files_R_environment_false_positive_blastncheck.Rdata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the R environments which contain all data used to create these summary results tables and figures. 

If you need to adjust the tables/figures you can easily load these documents into R and recreate the specific images you are after.



Summary_results_per_sample_contigs
----------------------------------

This folder contains sample specific results for the core module. A lot of this data is unnecessary to view as it has all been summarised in the folder 
**Summary_results_combined_figures_and_tables_contigs** but it also includes non summarised data including

**[Sample]_summarycontighits_assigned_assembly.txt** which contains all individual contig assignments and statistics around the BLAST assignments

as well as larger summary files like **[Sample]_top100Viralhits_contigs.txt** which contains the top 100 assignments instead of the top 10 reported in some of the combined summary top assignment tables

This folder also has a breakdown of the top returned SSU/LSU/CO1 reads hits per sample **[Sample]_top10returnedspecies_(SSU/LSU/CO1)**



raw_assemblies
--------------

This folder contains the generated assemblies for each sample. This can be useful if you want to do a deep dive, or you can take these contigs and put them into other
Analysis pipelines such as nf_core to allow for accurate bacterial identification and genome assembly without having to refilter and de-novo assemble raw data



CO1/LSU/SSU_results folders
---------------------------

When microbiome or host depletion modules are turned on three folders summarising the three marker genes are created. 


within each of these folders are three folders 

**blastn_contigs** which contains the BLASTn results of assembled contigs of these markers. 

**assembled_contigs** which contains the assembled contigs of each marker gene.

**read_homologies** which contains the results html figures of a Last common ancestor analysis used to identify what species are present.


.. image:: images/LCA_LSU_results.png
   :alt: MetaDIVE LSU_LCA_results_html
   :width: 300px
   :align: center





combined_raw_reads_and_contigs_viruses
--------------------------------------

This folder contains two key results types split by sample name. 


[Sample]virusall_sums.(csv/html)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files are very similar to the **Top_viral_hits_combined_raws_plus_contigs** having the same columns and data however it is sample specific results instead of combined.

This data is also in an html format for easier manual interrogation of important samples. Here the user can type any variables into the search bar
to look for patterns or key viral groups of interest.

.. image:: images/Results_table_html.png
   :alt: MetaDIVE Results table html
   :width: 300px
   :align: center


If the diverged reads/contigs module is run, a more complete csv file titled 

**[Sample]_virusall_sums_incl_diverged_reads.csv** whill be available. This file contains an additional 
column of data for additional diverged reads detected.


Individual viruses detected
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every virus identified in **[Sample]_virusall_sums** a folder is created for that specific virus. 
in this folder will be all contigs and all reads assigned to that viral species in fasta and fastq format for quick manual inspection
in programs like geneious. 



Metabat_binned_contigs
----------------------

To be completed




Reference_guided_viral_genomes
------------------------------

To be completed