library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(htmlwidgets)
library(magrittr)
library(argparse)
library(hrbrthemes)
library(htmltools)
library(hrbrthemes)
library(sankeyD3)
library(pavian)
library(phylotools)
library("taxonomizr")



if (!require("d3treeR", character.only = TRUE)) {
  devtools::install_github("timelyportfolio/d3treeR",upgrade="never",dependencies=FALSE)
  library("d3treeR")
} else {
  library("d3treeR")
}




parser <- ArgumentParser(description= 'Summarising results filtering and assembly')

parser$add_argument('--fastplog', '-a', help= 'Log for fastp')
parser$add_argument('--phiXlog', '-O', help= 'Log for Phix')
parser$add_argument('--CO1_bowtie', '-b', help= 'Log for CO1 bowtie alignment')
parser$add_argument('--LSU_bowtie', '-c', help= 'Log for LSU bowtie alignment')
parser$add_argument('--SSU_bowtie', '-d', help= 'Log for SSU bowtie alignment')
parser$add_argument('--CO1microbiome', '-e', help= 'results for CO1 LCA analysis')
parser$add_argument('--LSUmicrobiome', '-f', help= 'results for LSU LCA analysis')
parser$add_argument('--SSUmicrobiome', '-g', help= 'results for SSU LCA analysis')
parser$add_argument('--Assembly_log', '-i', help= 'Log for Assembly bowtie alignment')
parser$add_argument('--hostsp', '-j', help= ' Single top host species detected from LCA analysis of barcode markers')
parser$add_argument('--hostaligngenome', '-k', help= 'log for alignment of reads which aligned to downloaded host genome')
parser$add_argument('--diamondrawskingdoms', '-l', help= 'Summary of superkingdoms with all assigned raw reads')
parser$add_argument('--diamondrawslog', '-m', help= 'Diamond log for total assigned reads')
parser$add_argument('--krakenraws', '-n', help= 'Contigs file (kraken format) for raw reads assigned through Kraken2')
parser$add_argument('--rawtocontigsassignlog', '-o', help= 'Log for raw reads assigned back to generated contigs')
parser$add_argument('--samplename', '-p', help= 'sample name for given file/individual/pool')
parser$add_argument('--dohostdetect', '-q', help= 'Parameter, was host detect (and depletion) run')
parser$add_argument('--dodiamondraws', '-r', help= 'Parameter, was diamond run on raw reads')
parser$add_argument('--dokrakenraws', '-s', help= 'Parameter, was kraken2 run on raw reads')
parser$add_argument('--doblastn', '-t', help= 'Parameter, was Blastn run on contigs')
parser$add_argument('--domicrobiome', '-u', help= 'Parameter, were the microbiomes identified through rRNA and CO1')
parser$add_argument('--outputpath', '-v', help= 'Output path of the summary results')
parser$add_argument('--Log', '-w', help= 'log file to save results')
parser$add_argument('--programdir', '-x', help= 'base program directory')
parser$add_argument('--rawreadssamcontigs', '-y', help= 'sam file for raw reads assigned to contigs')
parser$add_argument('--Diamondtab', '-z', help= 'diamond contigs results table')
parser$add_argument('--Blastntab', '-A', help= 'Blastn contigs results table')
parser$add_argument('--rawssam', '-B', help= 'Raw sam results table')
parser$add_argument('--hostremovedcontigsfa', '-C', help= 'Proportion of contigs assigned to host')
parser$add_argument('--hostremovedcontigssam', '-D', help= 'Number of raw reads assigning to the host contigs filtered out')
parser$add_argument('--hostalignmentreads', '-E', help= 'Number of raw reads assigning to the host contigs in initial removal')
parser$add_argument('--dohostcontigsrawsalign', '-F', help= 'Number of raw reads assigning to contigs which aligned to the host genome')
parser$add_argument('--identifiedhost', '-G', help= 'The identified most likely host')
parser$add_argument('--Assemblyused', '-H', help= 'The assembly program used')
parser$add_argument('--hostgenomefailed', '-I', help= 'Was a host genome successfully found yes or no')
parser$add_argument('--Accnode', '-K', help= 'Accessiontaxa name_node filepath')
parser$add_argument('--dofalseposcontigschoice', '-L', help= 'Was a secondary blastn step taken to confirm viral origin')
parser$add_argument('--falseposcontigsrem', '-M', help= 'Viral contigs confirmed from blastn')
parser$add_argument('--doblastnassignmentsonly', '-N', help= 'whether to rename viral assignments from protein to nucleotide similarity')

xargs<- parser$parse_args()

# define base parameters and parameter variables for which analyses were run
NAMES <- xargs$samplename
basepath <- xargs$programdir
resultspath <- xargs$outputpath
outtablespath <- paste0(basepath,resultspath)
dohostdetect <- xargs$dohostdetect
dodiamondraws <- xargs$dodiamondraws
dokrakenraws <- xargs$dokrakenraws
dodoblastn <- xargs$doblastn
domicrobiome <- xargs$domicrobiome
docontigfalsepos <- xargs$dofalseposcontigschoice
docontigfalseposrenamespecies <- xargs$doblastnassignmentsonly
Assemblyused <- xargs$Assemblyused
Hostgenomefound <- xargs$hostgenomefailed


AccessionNamenode <- xargs$Accnode


# Starting analysis. Raw data filtering
# Each step here uses log files from bowtie alignments or other log files. 
# Each sample will have a log file with some values which are formatted consistently. As all steps are main steps, they can't be turned off
# and so if else or presence checks for each file are not needed
readLines(xargs$fastplog) -> fastpfile
length(fastpfile)
readLines(xargs$phiXlog) -> phixlogfile
length(phixlogfile)
readLines(xargs$CO1_bowtie)  -> bowtieCO1file
paste0(NAMES,"length bowtie CO1 ", length(bowtieCO1file))
readLines(xargs$LSU_bowtie) -> bowtieLSUfile
paste0(NAMES,"length bowtie LSU ", length(bowtieLSUfile))
readLines(xargs$SSU_bowtie) -> bowtieSSUfile
paste0(NAMES,"length bowtie SSU ", length(bowtieSSUfile))


# Function for reading in blast results tables which prevents errors associated with special characters in the species name
	replace_newline_within_quotes <- function(line) {
	  parts <- strsplit(line, "\"")[[1]]
	  for (i in seq_along(parts)) {
	    if (i %% 2 == 0) {
	      parts[[i]] <- gsub("\n", " ", parts[[i]])
	    }
	  }
	  paste(parts, collapse = "\"")
	}





summary_reads_table <- as.data.frame(matrix(nrow=1, ncol=12))

colnames(summary_reads_table) <- c("raw_reads","filtered_reads", "percentage_rem_QC", "Phix_filtered", "CO1_filtered", "Percentge_CO1","LSU_filtered", "Percentge_LSU","SSU_filtered", "Percentge_SSU", "remaining_reads","percentage_reads_remaining")

# Extract total count for raw reads
read1_line_index <- grep("Read1 before filtering", fastpfile)

# Extract the next line (which contains the total reads)
total_reads_line <- fastpfile[read1_line_index + 1]

raw1 <- as.numeric(sub(".*total reads: ([0-9]+).*", "\\1", total_reads_line))

read1_line_index <- grep("Read2 before filtering", fastpfile)

# Extract the next line (which contains the total reads)
total_reads_line <- fastpfile[read1_line_index + 1]

raw2 <- as.numeric(sub(".*total reads: ([0-9]+).*", "\\1", total_reads_line))


# Extract filtered reads number

read1_line_index <- grep("Read1 after filtering", fastpfile)

# Extract the next line (which contains the total reads)
total_reads_line <- fastpfile[read1_line_index + 1]

filtered1 <- as.numeric(sub(".*total reads: ([0-9]+).*", "\\1", total_reads_line))

read1_line_index <- grep("Read2 after filtering", fastpfile)

# Extract the next line (which contains the total reads)
total_reads_line <- fastpfile[read1_line_index + 1]

filtered2 <- as.numeric(sub(".*total reads: ([0-9]+).*", "\\1", total_reads_line))


total_raw <- ((raw1 + raw2))

# repeat for filtered data 
total_filtered <- ((filtered1 + filtered2))
total_filteredremoved <- (total_raw - total_filtered)

# calculate the percentage of read pairs which remain after filtering 
percentremain <- ((total_filtered / total_raw)*100)

summary_reads_table$raw_reads <- total_raw
summary_reads_table$filtered_reads <- total_filteredremoved
summary_reads_table$percentage_rem_QC  <- percentremain 

# 3. CO1 hits removed. Also adding in the Phix removal step here based on loss from prior fastp and total inputted into CO1 (can actually read on Phix if needed but its usually <10000 reads (<0.1%)
# and such a small component)



reads_line_index <- grep("reads; of these:", phixlogfile)

# Extract the line containing the reads information
reads_line <- phixlogfile[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_readpairs_prePhixfilter <- as.numeric(sub("([0-9]+) reads;.*", "\\1", reads_line))
total_reads_prePhixfilter <- (total_readpairs_prePhixfilter *2)



reads_line_index <- grep("reads; of these:", bowtieCO1file)

# Extract the line containing the reads information
reads_line <- bowtieCO1file[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_readpairs_preCO1filter <- as.numeric(sub("([0-9]+) reads;.*", "\\1", reads_line))
total_reads_preCO1filter <- (total_readpairs_preCO1filter *2)



reads_line_index <- grep("reads; of these:", bowtieLSUfile)
# Extract the line containing the reads information
reads_line <- bowtieLSUfile[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_readpairs_preLSUfilter <- as.numeric(sub("([0-9]+) reads;.*", "\\1", reads_line))
total_reads_preLSUfilter <- (total_readpairs_preLSUfilter *2)

# SSU is a bit weird because I have thrown the unassigned single reads back into the filter so they can be grabbed by a later raw reads check if required. 
# Ths means the total reads is not = 2* reads in first line of the log. 


reads_line_index <- grep(" were paired; of these", bowtieSSUfile)
# Extract the line containing the reads information
reads_line <- bowtieSSUfile[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_readpairs_preSSUfilter1 <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", reads_line))
total_reads_preSSUfilter1a <- (total_readpairs_preSSUfilter1 *2)

reads_line_index <- grep(" were unpaired", bowtieSSUfile)
# Extract the line containing the reads information
reads_line <- bowtieSSUfile[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_reads_preSSUfilter2 <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", reads_line))
total_reads_preSSUfilter <- (total_reads_preSSUfilter1a + total_reads_preSSUfilter2)



mates_line_index <- grep("mates make up the pairs", bowtieSSUfile)

# The target line is the one immediately after this line
target_line <- bowtieSSUfile[mates_line_index + 1]

# Use a regular expression to extract the first number in the target line
readsunalignedSSU1a <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", target_line))


mates_line_index <- grep(" were unpaired; of these", bowtieSSUfile)

# The target line is the one immediately after this line
target_line <- bowtieSSUfile[mates_line_index + 1]

# Use a regular expression to extract the first number in the target line
readsunalignedSSU1b <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", target_line))

readsunalignedSSU <- (readsunalignedSSU1a +readsunalignedSSU1b)


readsalignedSSU <- (total_reads_preSSUfilter - readsunalignedSSU )



alignment_rate_index <- grep("overall alignment rate", bowtieSSUfile)

# Extract the line containing the alignment rate information
alignment_rate_line <- bowtieSSUfile[alignment_rate_index]

# Use a regular expression to extract the percentage value
alignment_rate <- as.numeric(sub("([0-9\\.]+)% overall alignment rate.*", "\\1", alignment_rate_line))


readsfilteredPhix <- (total_reads_prePhixfilter-total_reads_preCO1filter)
readsfilteredCO1 <- (total_reads_preCO1filter-total_reads_preLSUfilter)
readsfilteredLSU <- (total_reads_preLSUfilter-total_reads_preSSUfilter)
readsfilteredSSU <- (readsalignedSSU)
readsfilteredSSU  <- round(readsfilteredSSU)
Percentage_CO1 <- ((readsfilteredCO1 / total_reads_preCO1filter)*100)
Percentage_LSU <- ((readsfilteredCO1 / total_reads_preCO1filter)*100)
Percentage_SSU <- alignment_rate

remainingreads <- readsunalignedSSU


summary_reads_table$Phix_filtered <- readsfilteredPhix 

summary_reads_table$CO1_filtered <- readsfilteredCO1 

summary_reads_table$Percentge_CO1 <- Percentage_CO1

summary_reads_table$LSU_filtered <- readsfilteredLSU

summary_reads_table$Percentge_LSU <- Percentage_LSU

summary_reads_table$SSU_filtered <- readsfilteredSSU

summary_reads_table$Percentge_SSU <- Percentage_SSU

summary_reads_table$remaining_reads <- remainingreads

summary_reads_table$percentage_reads_remaining <- ((remainingreads / summary_reads_table$raw_reads ) *100)

summary_reads_table$SAMPLE <- NAMES
write.table(summary_reads_table,file=(paste0(outtablespath,NAMES,"_summary_reads_filtering.txt")),sep="\t",row.names=FALSE)




#Now finished raw reads analysis 

# Assembly of contigs/reads general pipeline

if(dohostdetect =="yes") {
  
  if(Hostgenomefound == "yes") {
    # read lines for the sample 
    hostgenomecontigpresencehits <- length(readLines(xargs$hostremovedcontigsfa))
    if(hostgenomecontigpresencehits >=2) {
      hostassignedcontigsfa <- length(grep(">",xargs$hostremovedcontigsfa))
      numberhostcontigs <- hostassignedcontigsfa 
    }
    if(hostgenomecontigpresencehits < 1) {
      numberhostcontigs <- 0
    }
  }
  
  
  
  if(Hostgenomefound != "yes") {
    numberhostcontigs <- 0
  }
  
  
}





# Repeat for if no host filtering was done. Contig numbers then will be 0

if(dohostdetect =="no") {
  
  numberhostcontigs <- 0
  
}


rm(bowtieCO1file,bowtieLSUfile,bowtieSSUfile,fastpfile)

Diamondpresence <-readLines(xargs$Diamondtab)
if (length(Diamondpresence) >=1) {

	modified_lines <- lapply(Diamondpresence, replace_newline_within_quotes)
	modified_lines  <- unlist(modified_lines)
	modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
	# Read the modified lines using read.table
	Diamondhits <- read.table(text = modified_lines, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

	# Replace the unique character sequence back to '\n' within the dataframe
	Diamondhits[] <- lapply(Diamondhits, function(col) gsub("###NEWLINE###", "\n", col))

	rm(modified_lines)

#Diamondhits <- read.table(file=xargs$Diamondtab, sep="\t",header=TRUE,row.names=NULL, fill=TRUE,quote="")

  paste0(NAMES," dim Diamondhits ", dim(Diamondhits))

colnames(Diamondhits) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp", "multiplesp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")

Diamondhits$pident <- as.numeric(Diamondhits$pident)
Diamondhits$length<- as.numeric(Diamondhits$length)
Diamondhits$qcovhsp<- as.numeric(Diamondhits$qcovhsp)


# Diamond returns values as aa smiliarity, convert alignment lengths to nt
Diamondhits<- Diamondhits%>% 
  mutate(length = length * 3)



}
if (length(Diamondpresence) ==0) {
  Diamondhits <- as.data.frame(matrix(nrow=0,ncol=19))
  paste0(NAMES," Diamondx returned no findings")
}

blastnpresence <-readLines(xargs$Blastntab)
if (length(blastnpresence) >=1) {

	modified_lines <- lapply(blastnpresence, replace_newline_within_quotes)
	modified_lines  <- unlist(modified_lines)
	modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
	# Read the modified lines using read.table
	Blastnhits <- read.table(text = modified_lines, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

	# Replace the unique character sequence back to '\n' within the dataframe
	Blastnhits[] <- lapply(Blastnhits, function(col) gsub("###NEWLINE###", "\n", col))

	rm(modified_lines)


  #Blastnhits <- read.table(file=xargs$Blastntab, sep="\t",header=TRUE,row.names=NULL, fill=TRUE,quote="")
  paste0(NAMES," dim Blastnhits ", dim(Blastnhits))
}
if (length(blastnpresence) ==0) {
  Blastnhits <- as.data.frame(matrix(nrow=0,ncol=19))
  paste0(NAMES," Blastnhits returned no findings")
}

colnames(Blastnhits) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp", "multiplesp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")

Blastnhits$pident <- as.numeric(Blastnhits$pident)
Blastnhits$length<- as.numeric(Blastnhits$length)
Blastnhits$qcovhsp<- as.numeric(Blastnhits$qcovhsp)



# adding in check for if any viral contigs were returned after false positive check. This will get it to point where it can be read in to summarise all reads and raws. (Combined_assigned_contigsfp)

if ( docontigfalsepos == 'yes') {
  Diamondpresencefp <-readLines(xargs$falseposcontigsrem)
  if (length(Diamondpresencefp) >=1) {

	modified_lines <- lapply(Diamondpresencefp, replace_newline_within_quotes)
	modified_lines  <- unlist(modified_lines)
	modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
	# Read the modified lines using read.table
	Diamondhitsfp <- read.table(text = modified_lines, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

	# Replace the unique character sequence back to '\n' within the dataframe
	Diamondhitsfp[] <- lapply(Diamondhitsfp, function(col) gsub("###NEWLINE###", "\n", col))

	rm(modified_lines)


    #Diamondhitsfp <- read.table(file=xargs$falseposcontigsrem, sep="\t",header=TRUE,row.names=NULL, fill=TRUE,quote="")
    paste0(NAMES," dim Diamondhits after false positive check ", dim(Diamondhitsfp))
    
    colnames(Diamondhitsfp) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp", "multiplesp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies")
	Diamondhitsfp$pident <- as.numeric(Diamondhitsfp$pident)
	Diamondhitsfp$length<- as.numeric(Diamondhitsfp$length)
	Diamondhitsfp$qcovhsp<- as.numeric(Diamondhitsfp$qcovhsp)
   

	 Diamondhitsfp <- Diamondhitsfp%>% 
      		mutate(length = length * 3)
    
    
  }
  if (length(Diamondpresencefp) ==0) {
    Diamondhitsfp <- as.data.frame(matrix(nrow=0,ncol=19))
    paste0(NAMES," Diamondx after false positive check returned no findings")
    colnames(Diamondhitsfp) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp", "multiplesp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")
  }
  

  
  Combined_assigned_contigsfp <- rbind(Diamondhitsfp,Blastnhits)
  
}

Combined_assigned_contigs <- rbind(Diamondhits,Blastnhits)
Combined_assigned_contigspreblastnupdates <- rbind(Diamondhits,Blastnhits)
i=1
if (docontigfalseposrenamespecies == 'confirmed') {

	for (i in c(1:nrow(Diamondhitsfp))) {
		grep(paste0("^",Diamondhitsfp$qseqid[i],"$"),Combined_assigned_contigs$qseqid) -> idxval2
		
		if(length(idxval2) >=1) {

		Combined_assigned_contigs[idxval2,] <- Diamondhitsfp[i,]
		}
		# Adding new double check for if no idx values returned. This occurs only when
		# blastx identified nothing and kraken identified virus and blastn confirmed virus. Very uncommon but needed!
		# just rbind the new row in
		if(length(idxval2) ==0) {
			Combined_assigned_contigs <- rbind(Combined_assigned_contigs,Diamondhitsfp[i,])

		}

	}

}

summary_contigs_table <- as.data.frame(matrix(nrow=1, ncol=15))

colnames(summary_contigs_table) <- c("Number_contigs","Maximum_sized_contig", "N50_contigs", "Contigs_assigned_to_host", "Assigned_contigs_Diamond", "Assigned_contigs_Blastn","Percentage_of_contigs_assigned_total", "raw_reads_assigned_to_classified_contigs","Non_host_reads_assigned_via_protein_search","Non_host_reads_assigned_via_nucleotide_search", "Raw_reads_assigned_to_host_genome", "Raw_reads_assigned_to_host_aligned_contigs","Raw_reads_assigned_to_host_blastn_diamondx", "Raw_reads_assigned_raw_diamond", "raw_reads_not_assigned")

# Trinity doesn't produce a log file for assembly like Megahit does. Need to separate them here and gnerate the sum stats differently 
if ( Assemblyused =='Megahit') {
  
  readLines(xargs$Assembly_log)  -> contigsfile
  paste0(NAMES,"length contigslog ", length(contigsfile))
  
  # Will be parser sample name
  
  # 1. Contigs generated, min, max, n50
  
  # I am 99% sure that the rows are consistent in the log file and line 61 is always (assuming same Kmer ranges are used) the right row. I have written a grep check though to make sure
  # The phrase contigs is only used once in this file. (note, will need to confirm for trinity once this is run to
  grep("contigs,",contigsfile) ->idxval
  
  
  (str_extract(contigsfile[idxval], "[0-9]*.contigs")) -> contignum
  as.numeric(gsub(pattern= " contigs",replacement="", x=contignum)) -> contignum
  
  (str_extract(contigsfile[idxval], "max [0-9]*")) -> maxcontig
  as.numeric(gsub(pattern= "max ",replacement="", x=maxcontig)) -> maxcontig
  
  (str_extract(contigsfile[idxval], "N50 [0-9]*")) -> N50contig
  as.numeric(gsub(pattern= "N50 ",replacement="", x=N50contig)) -> N50contig
  
  summary_contigs_table$Number_contigs <- contignum
  summary_contigs_table$Maximum_sized_contig <- maxcontig
  summary_contigs_table$N50_contigs <- N50contig
  
}


if ( Assemblyused =='Trinity') {
  
  
  fasta <- read.fasta(xargs$Assembly_log)
  
  fasta$length <- NA
  
  for ( i in c(1:nrow(fasta))) {
    
    fasta$length[i] <- nchar(fasta$seq.text[i])
    
  }
  contignum <- nrow(fasta)
  summary_contigs_table$Number_contigs <- contignum
  summary_contigs_table$Maximum_sized_contig <- max(fasta$length)
  summary_contigs_table$N50_contigs <- median(fasta$length)
  
}




rm(Diamondpresence,blastnpresence)

summary_contigs_table$Contigs_assigned_to_host <- numberhostcontigs



# 2. Contigs assigned Diamond. 
# 3. Contigs assigned Blast. 
# 4. Contigs assigned total. 


readLines(xargs$rawtocontigsassignlog) -> bowtieraws
paste0(NAMES,"length rawbowtie ", length(bowtieraws))



rawssam <- read.table(xargs$rawssam, sep = "\t", fill = TRUE, row.names = NULL,header = FALSE)

paste0(NAMES,"dim rawssam ", dim(rawssam))

rawssam <- rawssam[,c(1,3)]

colnames(rawssam) <- c("read", "contig")

plyr::count(df = rawssam ,vars = "contig") -> freqsummary


freqsummary$contigassignment <- NA
freqsummary$superkingdom <- NA
freqsummary$phylum <- NA
freqsummary$class<- NA
freqsummary$order<- NA
freqsummary$family<- NA
freqsummary$genus<- NA
freqsummary$species<- NA
freqsummary$subspecies<- NA
freqsummary$percentident <- NA
freqsummary$contigalignlength <- NA

for (i in c(1:nrow(freqsummary))) {
  
  if(length(grep(paste0(freqsummary$contig[i],"$"),Diamondhits$qseqid))==1) {
    
    grep(paste0(freqsummary$contig[i],"$"), Diamondhits$qseqid) -> idxval
    
    freqsummary[i,3] <- "Diamondx"
    freqsummary[i,12] <- Diamondhits[idxval,3]
    freqsummary[i,13] <- Diamondhits[idxval,4]
    freqsummary[i,4:11] <- Diamondhits[idxval,12:19]
    
    
  }
  
  if(length(grep(paste0(freqsummary$contig[i],"$"),Blastnhits$qseqid))==1) {
    
    grep(paste0(freqsummary$contig[i],"$"), Blastnhits$qseqid) -> idxval
    
    freqsummary[i,3] <- "Blastn"
    freqsummary[i,12] <- Blastnhits[idxval,3]
    freqsummary[i,13] <- Blastnhits[idxval,4]
    freqsummary[i,4:11] <- Blastnhits[idxval,12:19]
    
  }
  
  
} 

cat(paste0("freqsummary worked", "\t"))
#Save the freqsum file of the diamond hits pre editing for blastn false positive (if option is chosen)
allassignedfreqspreblastnfpcheck <- subset(freqsummary,!is.na(freqsummary$contigassignment))
allassignedfreqspreblastnfpcheckNas <- freqsummary


i=1
if (docontigfalseposrenamespecies == 'confirmed') {
  if (nrow(Diamondhitsfp) >= 1) {
    for (i in c(1:nrow(Diamondhitsfp))) {
      
      idxval2 <- grep(paste0("^", Diamondhitsfp$qseqid[i], "$"), freqsummary$contig)
      
      # Check if idxval2 is not empty before accessing elements in freqsummary
      if (length(idxval2) > 0) {
        superkingombefore <- freqsummary$superkingdom[idxval2]
        superkingomafter <- Diamondhitsfp[i, 12]
        
        # Check if both superkingombefore and superkingomafter are non-empty and not NA
        if (!is.na(superkingombefore) && !is.na(superkingomafter) && superkingombefore != superkingomafter) {
          # Both values are not NA and they are different
          print(paste0("The superkingdoms of the blastx and blastn do not match."))
          print(paste0("Changing best blastx to blastn, resulting in ", freqsummary[idxval2, 11], " changing to ", Diamondhitsfp[i, 19]))
        }
        
        # Update freqsummary only if idxval2 is valid
        freqsummary[idxval2, 4:11] <- Diamondhitsfp[i, 12:19]
        freqsummary[idxval2, 12:13] <- Diamondhitsfp[i, 3:4]
        freqsummary[idxval2, 3] <- "Blastn"
      }
    }
  }
}






blastnassignedfreqs <- subset(freqsummary,freqsummary$contigassignment=="Blastn")
Diamondxassignedfreqs <- subset(freqsummary,freqsummary$contigassignment=="Diamondx")

# Adding in a section to get out only virus contigs to filter the raws sam before removing to extract all known virus reads for the contigs to feed into final viral section.

allassignedfreqs <- subset(freqsummary,!is.na(freqsummary$contigassignment))
allassignedfreqsvirus <- subset(allassignedfreqs,allassignedfreqs$superkingdom=="Viruses")
reads_contigs_virus <- rawssam[rawssam$contig %in% allassignedfreqsvirus$contig, ]

rm(rawssam)


if (dohostdetect =="yes") { 
  hostspecies <- read.table(file = xargs$identifiedhost, sep="\t",header=FALSE,row.names=NULL, fill=TRUE)
  
  hostspecies <- hostspecies[1,1]
  
  if (!(is.na(hostspecies))) {
    
    rows <- grep(hostspecies,allassignedfreqs$species)
    if (length(rows) >=1) {
      allassignedfreqshost <- allassignedfreqs[rows,]
      readsassignedclassifiedcontigshost <- sum(allassignedfreqshost$freq)
      ############
      allassignedfreqsnohost<- allassignedfreqs[-rows,]
      
      rows <- grep(hostspecies,blastnassignedfreqs$species)
      blastnfreqshost <- blastnassignedfreqs[rows,]
      readsassignedblastnhost <- sum(blastnfreqshost$freq)
      
      
      rows <- grep(hostspecies,Diamondxassignedfreqs$species)
      Diamondxfreqshost <- Diamondxassignedfreqs[rows,]
      readsassigneddiamondxhost <- sum(Diamondxfreqshost$freq)
    }
    if(length(rows) <1 ) {
      allassignedfreqshost <- allassignedfreqs[rows,]
      allassignedfreqsnohost <- allassignedfreqs
      
      
      blastnfreqshost <- blastnassignedfreqs[rows,]
      readsassignedblastnhost <- sum(blastnfreqshost$freq)
      Diamondxfreqshost <- Diamondxassignedfreqs[rows,]
      readsassigneddiamondxhost <- sum(Diamondxfreqshost$freq)
      
    }
    
    
  }
  
  if ((is.na(hostspecies))) {
    readsassignedclassifiedcontigshost <- 0
    readsassignedblastnhost <- 0
    readsassigneddiamondxhost <- 0
    allassignedfreqsnohost <- allassignedfreqs
  }
}
if (dohostdetect =="no") {
  
  readsassignedclassifiedcontigshost <- 0
  readsassignedblastnhost <- 0
  readsassigneddiamondxhost <- 0
  allassignedfreqsnohost <- allassignedfreqs
}


readsassignedclassifiedcontigs <- sum(allassignedfreqs$freq)
#readsassignedclassifiedcontigs <- (readsassignedclassifiedcontigs - readsassignedclassifiedcontigshost)

readsassignedblastn <- sum(blastnassignedfreqs$freq)
readsassignedblastn <- (readsassignedblastn - readsassignedblastnhost )

readsassigneddiamondx <- sum(Diamondxassignedfreqs$freq)
readsassigneddiamondx <- (readsassigneddiamondx - readsassigneddiamondxhost )


cat(paste0("counting read assignments to diamond and blast now completed", "\t"))

Diamondhitslen <- nrow(Diamondhits)
Blastnhitslen <- nrow(Blastnhits)

Percentagecontigsassigned <- ((( Diamondhitslen + Blastnhitslen + numberhostcontigs) / contignum ) *100) 


summary_contigs_table$Assigned_contigs_Diamond <- Diamondhitslen 
summary_contigs_table$Assigned_contigs_Blastn <- Blastnhitslen
summary_contigs_table$Percentage_of_contigs_assigned_total <- Percentagecontigsassigned



length(bowtieraws) -> X
percentage <- as.numeric(str_extract(bowtieraws[X], "[0-9]+\\.[0-9]+"))
percentage_unassigned <- (100-percentage)

reads_line_index <- grep("reads; of these:", bowtieraws)

# Extract the line containing the reads information
reads_line <- bowtieraws[reads_line_index]

# Use a regular expression to extract the number before "reads"
total_readpairs_prealigntoraws <- as.numeric(sub("([0-9]+) reads;.*", "\\1", reads_line))
total_reads_prealigningtoraws <- (total_readpairs_prealigntoraws *2)



unaligned <- (((percentage_unassigned / 100)*total_reads_prealigningtoraws ))

as.numeric(str_extract(bowtieraws[1], "[0-9]+")) -> total

aligned_raws <- ( total_reads_prealigningtoraws - unaligned )


# Not sure if I want the classified contigs one
summary_contigs_table$raw_reads_assigned_to_classified_contigs <- (readsassignedclassifiedcontigs)
summary_contigs_table$Non_host_reads_assigned_via_protein_search <- (readsassigneddiamondx)
summary_contigs_table$Non_host_reads_assigned_via_nucleotide_search <- (readsassignedblastn)
summary_contigs_table$Raw_reads_assigned_to_host_blastn_diamondx <- (readsassigneddiamondxhost + readsassignedblastnhost)
# Potentially rename to raw reads not assigned total, then have grey reads as measure here. 
#summary_contigs_table$raw_reads_not_assigned_to_contigs <- ( summary_reads_table$remaining_reads - aligned_raws )



if (dohostdetect =="yes") {
  if(Hostgenomefound == "yes") {
    # read lines for the sample 
    
    
    readLines(xargs$hostaligngenome) -> hostrawsalign
    
    readLines(xargs$dohostcontigsrawsalign) -> readstohostcontigsalign
    
    if (length(hostrawsalign) >=1) {




	reads_line_index <- grep(" were paired; of these", hostrawsalign)
	# Extract the line containing the reads information
	reads_line <- hostrawsalign[reads_line_index]

	# Use a regular expression to extract the number before "reads"
	total_readpairss_prehostfilter1 <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", reads_line))
	total_reads_prehostfilter1a <- (total_readpairss_prehostfilter1 *2)

	reads_line_index <- grep(" were unpaired", hostrawsalign)
	# Extract the line containing the reads information
	reads_line <- hostrawsalign[reads_line_index]

	# Use a regular expression to extract the number before "reads"
	total_readpairs_prehostfilter2 <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", reads_line))
	total_reads_prehostfilter <- (total_reads_prehostfilter1a + total_readpairs_prehostfilter2)



	mates_line_index <- grep("mates make up the pairs", hostrawsalign)

# The target line is the one immediately after this line
	target_line <- hostrawsalign[mates_line_index + 1]

#Use a regular expression to extract the first number in the target line
	readsunalignedhost1a <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", target_line))


	mates_line_index <- grep(" were unpaired; of these", hostrawsalign)
# The target line is the one immediately after this line
	target_line <- hostrawsalign[mates_line_index + 1]

# Use a regular expression to extract the first number in the target line
	readsunalignedhost1b <- as.numeric(sub("^\\s*([0-9]+) \\(.*", "\\1", target_line))

	readsunalignedhost <- (readsunalignedhost1a + readsunalignedhost1b)


	readsalignedhost <- (total_reads_prehostfilter - readsunalignedhost )






      total <- total_reads_prehostfilter
      length(hostrawsalign) -> X
      percentage <- as.numeric(str_extract(hostrawsalign[X], "[0-9]+\\.[0-9]+"))
      percentage_unassigned <- (100-percentage)
      unassignedrawhost <-  readsunalignedhost     
      # Put inside first because if the sample doesn't get any raw read hits to a host genome its nor
      # going to get contig hits later
      if (length(readstohostcontigsalign) >=1) {
        
        totalhostcont <- (as.numeric(str_extract(readstohostcontigsalign[1], "[0-9]+"))*2)
        length(readstohostcontigsalign) -> X
        percentage <- as.numeric(str_extract(readstohostcontigsalign[X], "[0-9]+\\.[0-9]+"))
        percentage_unassigned <- (100-percentage)
        unassignedhostcont <- ((percentage_unassigned*totalhostcont)/100)
        
      }
      
      if (length(readstohostcontigsalign) ==0) {
        
        totalhostcont <- unassignedrawhost
        length(hostrawsalign) -> X
        percentage <- as.numeric(str_extract(hostrawsalign[X], "[0-9]+\\.[0-9]+"))
        percentage_unassigned <- (100-percentage)
        unassignedhostcont <- unassignedrawhost
        
        unassignedhostcont  <- round(unassignedhostcont)
      }
      
      
    }
    
    if (length(hostrawsalign) ==0) {
      
      unassigned <- summary_reads_table$remaining_reads
      total <- summary_reads_table$remaining_reads
      
      unassignedhostcont <- summary_reads_table$remaining_reads
      unassignedrawhost<- summary_reads_table$remaining_reads
      totalhostcont <-  summary_reads_table$remaining_reads

	unassignedhostcont  <- round(unassignedhostcont)
	unassignedrawhost<- round(unassignedrawhost)
      
    }
    
    
    hostassignedreads <- (total - unassignedrawhost )
    hostassignedreads <- round(hostassignedreads)
    hostcontigsassignedreads <- (totalhostcont - unassignedhostcont )
    hostcontigsassignedreads <- round(hostcontigsassignedreads)
    if (length(readstohostcontigsalign) ==0){
      hostcontigsassignedreads <- 0
    }
    
    
    
  }
  
  if(Hostgenomefound == "no") {
    hostassignedreads <- 0
    hostcontigsassignedreads  <- 0
  }
  
  
}
# Repeat for if no host filtering was done. Contig numbers then will be 0

if (dohostdetect =="no") {
  
  hostassignedreads <- 0
  hostcontigsassignedreads  <- 0
}

summary_contigs_table$Raw_reads_assigned_to_host_genome <- hostassignedreads 
summary_contigs_table$Raw_reads_assigned_to_host_aligned_contigs <- hostcontigsassignedreads


cat(paste0("putting all summary table values together worked", "\t"))



# finished base contig identifications




#generate top hits and average identity for contigs
species_idvec <- unique(allassignedfreqsnohost$species)

contigsspecies <- as.data.frame(matrix(nrow=length(species_idvec), ncol=13))
colnames(contigsspecies) <- c("Frequency","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident","length")


cat(paste0("creating species idvec worked", "\t"))

freqsummaryonlyna <- subset(freqsummary, is.na(freqsummary$percentident))
readsfromunassignedcontigs <- sum(freqsummaryonlyna$freq)

freqsummarynona <- subset(freqsummary, !is.na(freqsummary$percentident))
freqsummarynona$percentident <- as.numeric(freqsummarynona$percentident)
for (i in c(1:length(species_idvec))) {
  
  contigssubset<- subset(freqsummarynona,freqsummarynona$species == species_idvec[i])
  
  contigsspecies[i,2:9] <- contigssubset[1,4:11]
  contigsspecies[i,1] <- sum(contigssubset$freq)
  
  if (!is.na(contigssubset$percentident[1])) {
    contigsspecies[i,10] <- mean(contigssubset$percentident)
    contigsspecies[i,11] <- min(contigssubset$percentident)
    contigsspecies[i,12] <- max(contigssubset$percentident)
    contigsspecies[i,13] <- mean(contigssubset$contigalignlength)
  }
}




cat(paste0("creating contig subsets worked", "\t"))




Bacteriaspecies <- subset(contigsspecies,contigsspecies$superkingdom=="Bacteria")
Viralspecies <- subset(contigsspecies,contigsspecies$superkingdom=="Viruses")
Eukaryotesspecies <- subset(contigsspecies,contigsspecies$superkingdom=="Eukaryota")

Bacteriasordered <- Bacteriaspecies[order(-Bacteriaspecies$Frequency),]

if (nrow(Bacteriasordered)<= 20) {
  
  Bacteriatop20 <- Bacteriasordered
  
}


if (nrow(Bacteriasordered) > 20) {
  
  Bacteriatop20 <- Bacteriasordered[1:20,]
  
}

if (nrow(Bacteriasordered) > 100) {
  
  Bacteriatop100 <- Bacteriasordered[1:100,]
  
}


if (nrow(Bacteriasordered)<= 100) {
  
  Bacteriatop100<- Bacteriasordered
  
}



Viralspeciesordered <- Viralspecies[order(-Viralspecies$Frequency),]

if (nrow(Viralspeciesordered)<= 20) {
  
  Viraltop20 <- Viralspeciesordered
  
}


if (nrow(Viralspeciesordered) > 20) {
  
  Viraltop20 <- Viralspeciesordered[1:20,]
  
}


if (nrow(Viralspeciesordered) > 100) {
  
  Viraltop100 <- Viralspeciesordered[1:100,]
  
}


if (nrow(Viralspeciesordered)<= 100) {
  
  Viraltop100 <- Viralspeciesordered
  
}



Eukaryotesspeciesordered <- Eukaryotesspecies[order(-Eukaryotesspecies$Frequency),]

if (nrow(Eukaryotesspeciesordered)<= 20) {
  
  Eukaryotestop20 <- Eukaryotesspeciesordered
  
}


if (nrow(Eukaryotesspeciesordered) > 20) {
  
  Eukaryotestop20 <- Eukaryotesspeciesordered[1:20,]
  
}


if (nrow(Eukaryotesspeciesordered) > 100) {
  
  Eukaryotestop100 <- Eukaryotesspeciesordered[1:100,]
  
}


if (nrow(Eukaryotesspeciesordered)<= 100) {
  
  Eukaryotestop100<- Eukaryotesspeciesordered
  
}



if (nrow(Bacteriatop20) >=1) {
  Bacteriatop20$Sample <- NAMES
}
if (nrow(Eukaryotestop20) >=1) {
  Eukaryotestop20$Sample <- NAMES
}

if (nrow(Viraltop20) >=1) {
  Viraltop20$Sample <- NAMES
}

if (nrow(freqsummary) >=1) {
  freqsummary$Sample <- NAMES
}

if (nrow(allassignedfreqs) >=1) {
  allassignedfreqs$Sample <- NAMES
}



# Save all files generated
# Note I am currently investigating, but it appears that Diamondx and blastn perform poorly when larger queries are given (over 10k bp) 
write.table(Bacteriatop20,file=(paste0(outtablespath,NAMES,"_top20bacterialhits_contigs.txt")),sep="\t",row.names=FALSE)
write.table(Eukaryotestop20,file=(paste0(outtablespath,NAMES,"_top20Eukaryotehits_contigs.txt")),sep="\t",row.names=FALSE)
write.table(Viraltop20,file=(paste0(outtablespath,NAMES,"_top20Viralhits_contigs.txt")),sep="\t",row.names=FALSE)
write.table(freqsummary,file=(paste0(outtablespath,NAMES,"_summarycontighits_assembly.txt")),sep="\t",row.names=FALSE)
write.table(allassignedfreqs,file=(paste0(outtablespath,NAMES,"_summarycontighits_assigned_assembly.txt")),sep="\t",row.names=FALSE)
write.table(Bacteriatop100,file=(paste0(outtablespath,NAMES,"_top100bacterialhits_contigs.txt")),sep="\t",row.names=FALSE)
write.table(Eukaryotestop100,file=(paste0(outtablespath,NAMES,"_top100Eukaryotehits_contigs.txt")),sep="\t",row.names=FALSE)
write.table(Viraltop100,file=(paste0(outtablespath,NAMES,"_top100Viralhits_contigs.txt")),sep="\t",row.names=FALSE)



# Generate output_taxids for results
# Note it is currently pulling top 20 despite the name top10


Viraltop10 <- Viraltop20[1:20,]

Viraltop10 <- subset(Viraltop10,Viraltop10$average_percent_ident >=40 & Viraltop10$Frequency >=50)


if (nrow(Viraltop10) ==0) {
  Viraltop10 <- Viraltop20[1:20,]
  
  Viraltop10 <- subset(Viraltop10,Viraltop10$average_percent_ident >=30 & Viraltop10$Frequency >=20)
  
  
}





Viraltoptaxids <- as.data.frame(matrix(nrow=nrow(Viraltop10),ncol=1))
if (nrow(Viraltop10) !=0) {
  
  for ( i in c(1:nrow(Viraltop10))) {
    
    grep(Viraltop10$species[i],Combined_assigned_contigs$species) -> idx
    value <- Combined_assigned_contigs[idx[1],]
    
    Viraltoptaxids[i,1] <- value$staxidreduced
    
  }
  
}



# I want to add the viruses that show >80% similarity as well as top 10. 

#first. remove top 10 which will always be present
# Next subset by pwident 
Viraltop_idents <- subset(Viraltop100,Viraltop100$average_percent_ident >=80 & Viraltop100$Frequency >=50)
Viraltop_idents <- Viraltop_idents[order(-Viraltop_idents$Frequency), ]

Viraltop_identstaxids <- as.data.frame(matrix(nrow=nrow(Viraltop_idents),ncol=1))



if ( nrow(Viraltop_identstaxids) >=1) {
  
  for ( i in c(1:nrow(Viraltop_idents))) {
    
    grep(Viraltop_idents$species[i],Combined_assigned_contigs$species) -> idx
    value <- Combined_assigned_contigs[idx[1],]
    
    Viraltop_identstaxids[i,1] <- value$staxidreduced
    
  }
  
  Viraltoptaxids <- rbind(Viraltoptaxids,Viraltop_identstaxids)
  Viraltoptaxids <- Viraltoptaxids[!duplicated(Viraltoptaxids$V1), ]
  Viraltoptaxids <- as.data.frame(Viraltoptaxids)
  colnames(Viraltoptaxids) <- "V1"
  
}


if (nrow(Viraltop10) >=1) {
  
  if (nrow(Viraltoptaxids) >=1) {
    
    if (!is.na(Viraltop10[1,1])) {
      
      if (nrow(Viraltoptaxids) <=20) {
        
        write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE, col.names=FALSE,quote=FALSE)
      }
      if (nrow(Viraltoptaxids ) >20) {
        Viraltoptaxids <- Viraltoptaxids[1:20,]
        write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE, col.names=FALSE,quote=FALSE)
      }
      
      
    }
  }
  
}

if (is.na(Viraltop10[1,1])) {
  Viraltoptaxids <- as.data.frame(matrix(nrow=1,ncol=1))
  
  write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

if (!exists("Viraltoptaxids")) {
  Viraltoptaxids <- as.data.frame(matrix(nrow=1,ncol=1))
  
  write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}


if (nrow(Viraltop10) ==0) {
  Viraltoptaxids <- as.data.frame(matrix(nrow=1,ncol=1))
  
  write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}





if (dohostdetect=="yes") {
  
  if (nrow(allassignedfreqsnohost) >=1) {
    allassignedfreqsnohost$Sample <- NAMES
  }
  
  
  
  write.table(allassignedfreqsnohost,file=(paste0(outtablespath,NAMES,"_summarycontighits_assigned_no_host_assembly.txt")),sep="\t",row.names=FALSE)
}



# identifiedhost


# Now analysing the raw reads alignments to Diamond/Kraken 





summary_raw_reads_assignments <- as.data.frame(matrix(nrow=1, ncol=2))

colnames (summary_raw_reads_assignments) <- c("Number_of_raw_reads_assigned_through_Diamond", "Number_of_raw_reads_assigned_through_Kraken")


if (dodiamondraws == 'yes') {
  
  Diamondrawpresence <-readLines(xargs$diamondrawskingdoms)
  if (length(Diamondrawpresence >=1)) {
    
	modified_lines <- lapply(Diamondrawpresence, replace_newline_within_quotes)
	modified_lines  <- unlist(modified_lines)
	modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
	# Read the modified lines using read.table
	Diamondrawshits <- read.table(text = modified_lines, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

	# Replace the unique character sequence back to '\n' within the dataframe
	Diamondrawshits[] <- lapply(Diamondrawshits, function(col) gsub("###NEWLINE###", "\n", col))

	rm(modified_lines)

    paste0(NAMES," dim Diamondrawshits ", dim(Diamondrawshits))
    
    
    colnames(Diamondrawshits) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies")
    
  }
  if (length(Diamondrawpresence) ==0) {
    Diamondrawshits <- as.data.frame(matrix(nrow=0,ncol=18))
    paste0(NAMES," Diamond raws analysis returned no findings")
  }
  
  colnames(Diamondrawshits) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies")
  
  # Diamond returns length alignments as aa length, convert to nucleotide length
  
Diamondrawshits$pident <- as.numeric(Diamondrawshits$pident)
Diamondrawshits$length<- as.numeric(Diamondrawshits$length)
Diamondrawshits$qcovhsp<- as.numeric(Diamondrawshits$qcovhsp)


  Diamondrawshits<- Diamondrawshits%>% 
    mutate(length = length * 3)
  

  
  # 
  #Diamondrawslog <- readLines(xargs$diamondrawslog)
  
  #grep("queries aligned", Diamondrawslog) -> idx
  
  #readstotal <- Diamondrawslog[idx]
  
  #readstotal <- gsub(pattern =  " queries aligned.", replacement = "" ,x=readstotal)
  
  #readstotal <- as.numeric(readstotal)
  #readsassigned <- sum(Diamondrawskingdoms$freq)
  
  #summary_raw_reads_assignments$Number_of_raw_reads_assigned_through_Diamond <- readsassigned 
  
  # So the raw reads analysis in Diamond is blasting individual reads. This means it is possible for 1 read to assign
  # but its pair not to. All other analyses count read pairs which messes with the final counts of how to incorporate these
  # I can 1. double the remaining reads after all removals instead of Diamond or I can subtract the uniques of Diamond 
  # under the assumption that if one of the pair was subtracted the other should probably have been subtracted to. 
  # Doing the latter here because I don't like having everything in paired ends except the very last stat

  length(unique(Diamondrawshits$qseqid)) -> rawdiamondassignedreads
  
  summary_contigs_table$Raw_reads_assigned_raw_diamond <- rawdiamondassignedreads
  
  paste0(" No species level information for raw Diamond blast will be returned for low identity reads here, for details of assignments see directory /11_Diamond_ASSIGNED_RAWS/abundances or read the final output tables combining raws and contigs")
  
  summary_contigs_table$raw_reads_not_assigned <- (summary_reads_table$remaining_reads-(summary_contigs_table$Non_host_reads_assigned_via_protein_search + summary_contigs_table$Non_host_reads_assigned_via_nucleotide_search + summary_contigs_table$Raw_reads_assigned_to_host_genome + summary_contigs_table$Raw_reads_assigned_to_host_aligned_contigs + summary_contigs_table$Raw_reads_assigned_to_host_blastn_diamondx + summary_contigs_table$Raw_reads_assigned_raw_diamond + readsfromunassignedcontigs))
  
  summary_contigs_table$raw_reads_not_assigned <- round(summary_contigs_table$raw_reads_not_assigned)
  
  Diamondrawshitshighacc <- subset(Diamondrawshits,Diamondrawshits$pident >=55 & Diamondrawshits$length>=40)
  
  Diamondrawshitshighaccfreqs <- Diamondrawshitshighacc %>%
    group_by(species) %>%
    summarize(count = n())
  

  
  Diamondrawshitshighaccfreqs<- subset(Diamondrawshitshighaccfreqs,!(is.na(Diamondrawshitshighaccfreqs$species)))
  
  Diamondrawshitshighaccfreqs$superkingdom <- NA 
  Diamondrawshitshighaccfreqs$family <- NA 
  Diamondrawshitshighaccfreqs$taxid <- NA
  Diamondrawshitshighaccfreqs$min_aligned_length <- NA 
  Diamondrawshitshighaccfreqs$max_aligned_length <- NA 
  Diamondrawshitshighaccfreqs$mean_aligned_length <- NA 
  Diamondrawshitshighaccfreqs$min_identity <- NA 
  Diamondrawshitshighaccfreqs$max_identity <- NA 
  Diamondrawshitshighaccfreqs$mean_identity <- NA 
  
  if(nrow(Diamondrawshitshighaccfreqs) >=1) {
    for (l in c(1:nrow(Diamondrawshitshighaccfreqs) )) {
      
      Diamondrawshitshighaccspeciessubset <- Diamondrawshitshighacc %>%
        filter(grepl(Diamondrawshitshighaccfreqs[l,1], Diamondrawshitshighacc$species))
      
      
      
      Diamondrawshitshighaccfreqs$superkingdom[l] <- Diamondrawshitshighaccspeciessubset$superkingdom[1]
      Diamondrawshitshighaccfreqs$family[l] <- Diamondrawshitshighaccspeciessubset$family[1]
      Diamondrawshitshighaccfreqs$taxid[l] <- Diamondrawshitshighaccspeciessubset$staxids[1]
      Diamondrawshitshighaccfreqs$min_aligned_length[l] <- min(Diamondrawshitshighaccspeciessubset$length)
      Diamondrawshitshighaccfreqs$max_aligned_length[l] <- max(Diamondrawshitshighaccspeciessubset$length)
      Diamondrawshitshighaccfreqs$mean_aligned_length[l] <- mean(Diamondrawshitshighaccspeciessubset$length)
      Diamondrawshitshighaccfreqs$min_identity[l] <- min(Diamondrawshitshighaccspeciessubset$pident)
      Diamondrawshitshighaccfreqs$max_identity[l] <- max(Diamondrawshitshighaccspeciessubset$pident)
      Diamondrawshitshighaccfreqs$mean_identity[l] <- mean(Diamondrawshitshighaccspeciessubset$pident)
      
    }
    
    # This should become a variable to play around with! I have set it to 1 right now as I am not sure about the benefit of not having it at all given there will be externam filters the user can apply on the final result
    # Maybe its worth having for speed up purposes but it might be leading to a drop in true assignment rates

    Diamondrawshitshighaccfreqssignificant <- subset(Diamondrawshitshighaccfreqs,Diamondrawshitshighaccfreqs$count>=1)
    
    # Now need to fill in the taxids and then split into raw files and add viruses on
    
    
    for (i in c(1:nrow(Diamondrawshitshighaccfreqssignificant))) {
      
      if (is.na(Diamondrawshitshighaccfreqssignificant$taxid[i])) {
        
        Diamondrawshitshighaccfreqssignificant$taxid[i] <- taxonomizr::getId(Diamondrawshitshighaccfreqssignificant$species[i], sqlFile=AccessionNamenode)
        
      }
      
      
    }

    DiamondrawshitshighaccfreqssignificantEukaryotes <- subset(Diamondrawshitshighaccfreqssignificant,Diamondrawshitshighaccfreqssignificant$superkingdom=="Eukaryota")
    DiamondrawshitshighaccfreqssignificantBacteria <- subset(Diamondrawshitshighaccfreqssignificant,Diamondrawshitshighaccfreqssignificant$superkingdom=="Bacteria")
    DiamondrawshitshighaccfreqssignificantVirus<- subset(Diamondrawshitshighaccfreqssignificant,Diamondrawshitshighaccfreqssignificant$superkingdom=="Viruses")
    
    
    write.table(DiamondrawshitshighaccfreqssignificantBacteria,file=(paste0(outtablespath,NAMES,"bacterialhits_raw_reads.txt")),sep="\t",row.names=FALSE)
    write.table(DiamondrawshitshighaccfreqssignificantEukaryotes,file=(paste0(outtablespath,NAMES,"Eukaryotehits_raw_reads.txt")),sep="\t",row.names=FALSE)
    write.table(DiamondrawshitshighaccfreqssignificantVirus,file=(paste0(outtablespath,NAMES,"Viralhits_raw_reads.txt")),sep="\t",row.names=FALSE)
    
    
  }

  
  if (nrow(Diamondrawshitshighaccfreqs) <1) {
    
    # Create empty files to have as outputs
    file_bact_empty <- paste0(outtablespath,NAMES,"bacterialhits_raw_reads.txt")
    file.create(file_bact_empty)
    file_euk_empty <- paste0(outtablespath,NAMES,"Eukaryotehits_raw_reads.txt")
    file.create(file_euk_empty)
    file_vir_empty <- paste0(outtablespath,NAMES,"Viralhits_raw_reads.txt")
    file.create(file_vir_empty )
    
  }

  if (nrow(Diamondrawshitshighaccfreqs) >=1) {
    Viralraw_reads <- Diamondrawshitshighacc[Diamondrawshitshighacc$species %in% DiamondrawshitshighaccfreqssignificantVirus$species, ]
    # Can put secondary more stringent filter here for reducing impact of raw reads assignment? 
    Viralraw_reads <- subset(Viralraw_reads,Viralraw_reads$length>60 & Viralraw_reads$pident>50)
    
    Viralraw_readsnames <- Viralraw_reads$qseqid
    
    
    write.table(Viralraw_reads,file=(paste0(outtablespath,NAMES,"_virus_hits_all_details_results_table.txt")),sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)
    
    
    write.table(Viralraw_readsnames,file=(paste0(outtablespath,NAMES,"_raw_read_names_to_virus.txt")),sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)
    

    
    DiamondrawshitshighaccfreqssignificantVirusforlocalallignments <- subset(DiamondrawshitshighaccfreqssignificantVirus,DiamondrawshitshighaccfreqssignificantVirus$mean_identity >=80)
    DiamondrawshitshighaccfreqssignificantVirusforlocalallignments <- subset(DiamondrawshitshighaccfreqssignificantVirus,DiamondrawshitshighaccfreqssignificantVirus$count>=80)
    
    
    
    if (is.matrix(Viraltoptaxids) || (is.data.frame(Viraltoptaxids) && ncol(Viraltoptaxids) > 1)) {
      # Subset if Viraltoptaxids is a multi-column object
      Viraltoptaxids <- subset(Viraltoptaxids, !(is.na(Viraltoptaxids[,1])))
    } 
    if (is.vector(Viraltoptaxids)) {
      Viraltoptaxids <- as.data.frame(Viraltoptaxids)
      Viraltoptaxids <- as.matrix(Viraltoptaxids,ncol=1)
      # Leave it as is if it's already a vector
      
    }
    
    viral_rawstaxids <- DiamondrawshitshighaccfreqssignificantVirusforlocalallignments$taxid
    viral_rawstaxids <- as.numeric(viral_rawstaxids)
    viral_rawstaxids <- as.matrix(viral_rawstaxids,ncol=1)
    
    Viraltoptaxids <- rbind(Viraltoptaxids,viral_rawstaxids)
    
    Viraltoptaxids  <- unique(Viraltoptaxids[, 1])
    
    
    
    
  }
  
}

if (nrow(Diamondrawshitshighaccfreqs) <1) {
  file_virreads_empty <- paste0(outtablespath,NAMES,"_virus_hits_all_details_results_table.txt")
  file.create(file_virreads_empty)
  file_virnames_empty <- paste0(outtablespath,NAMES,"_raw_read_names_to_virus.txt")
  file.create(file_virnames_empty)
  
  
  Viraltoptaxids <- subset(Viraltoptaxids,!(is.na(Viraltoptaxids[,1])))
  
}

if (dodiamondraws == 'no') {
  summary_contigs_table$raw_reads_not_assigned <- (summary_reads_table$remaining_reads-(summary_contigs_table$Non_host_reads_assigned_via_protein_search + summary_contigs_table$Non_host_reads_assigned_via_nucleotide_search + summary_contigs_table$Raw_reads_assigned_to_host_genome + summary_contigs_table$Raw_reads_assigned_to_host_aligned_contigs + summary_contigs_table$Raw_reads_assigned_to_host_blastn_diamondx + readsfromunassignedcontigs))
}





if (!is.na(Viraltop10[1,1])) {
  
  write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE, col.names=FALSE,quote=FALSE)
}

if (is.na(Viraltop10[1,1])) {
  Viraltoptaxids <- as.data.frame(matrix(nrow=1,ncol=1))
  
  Viraltoptaxids <- subset(Viraltoptaxids,!(is.na(Viraltoptaxids[,1])))
  
  write.table(Viraltoptaxids,file=(paste0(outtablespath,NAMES,"_top10_Taxids_Viral_contigs.txt")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}






summary_contigs_table$SAMPLE <- NAMES
write.table(summary_contigs_table,file=(paste0(outtablespath,NAMES,"_summarycontigs_assembly_values.txt")),sep="\t",row.names=FALSE)

##############Up to here for building redundancy. Need to check whether Diamondtab is correct here. 

if (dokrakenraws == 'yes') {
  
  krakenraws <- pavian::read_report(xargs$krakenraws)
  
  
  if (!(is.null(dokrakenraws))) {
    
    
    virusonlyidx <- grep("d_Viruses",krakenraws$taxLineage)
    virusonlykraken <- krakenraws[virusonlyidx,]
    
    virusonlykrakenreads <- sum(virusonlykraken$taxonReads)
    summary_raw_reads_assignments$Number_of_raw_reads_assigned_through_Kraken <- virusonlykrakenreads
    
    
  }
  
  
  if ((is.null(dokrakenraws))) {
    
    paste0("Kraken returned zero assignments to species for any virus. Report cannot be generated")
    
    summary_raw_reads_assignments$Number_of_raw_reads_assigned_through_Kraken <- 0
    
  }
  
  
}


if (dodiamondraws == 'yes') {
  
  summary_raw_reads_assignments$Number_of_raw_reads_assigned_through_Diamond <- rawdiamondassignedreads
  
  if (nrow(Diamondrawshits) >=1) {
    Diamondrawshits$Sample <- NAMES
  }
  
  
  if (dohostdetect =="yes") {
    
    if( hostspecies[1] !="NA" ) { 
      
      rows <- grep(hostspecies,Diamondrawshits$species)
      if ( length(rows) >=1) {
        Diamondrawshits<- Diamondrawshits[-rows,]
      }
    }
  }
  
  if (dohostdetect =="no") {
    
  }
  
  
  write.table(Diamondrawshits,file=(paste0(outtablespath,NAMES,"_raw_diamond_hits.txt")),sep="\t",row.names=FALSE)
  
  
}

if (dokrakenraws == 'no') {
  
  summary_raw_reads_assignments$Number_of_raw_reads_assigned_through_Kraken <- 0
  
}

# Specialised pathways: Host identification and Microbiomes 

# CO1

CO1microbiome <- pavian::read_report(xargs$CO1microbiome)

if (!(is.null(CO1microbiome))) {
  
  topspecies <- subset(CO1microbiome, CO1microbiome$taxRank=="S")
  
  topspeciesordered <- topspecies[order(-topspecies$taxonReads),]
  
  if (nrow(topspeciesordered) >=1) {
    
    if (nrow(topspeciesordered) > 10) {
      topspeciestop10 <- topspeciesordered[1:10,]
    }
    if (nrow(topspeciesordered) <= 10) {
      topspeciestop10 <- topspeciesordered
    }
    
    topCO1hits <- as.data.frame(matrix(nrow=nrow(topspeciestop10), ncol=5))
    colnames(topCO1hits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_CO1_reads")
    
    strings <- str_extract(topspeciestop10$taxLineage, "f_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="f_",replacement="")
    topCO1hits$Family <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "g_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="g_",replacement="")
    topCO1hits$Genus <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "s_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="s_",replacement="")
    topCO1hits$Species <- unlist(strings)
    
    topCO1hits$read_hits <- topspeciestop10$taxonreads
    topCO1hits$percentage_of_CO1_reads <- topspeciestop10$percentage
    
  }
  
  if (nrow(topspecies) ==0) { 
    
    topCO1hits <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(topCO1hits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_CO1_reads")
    
    topCO1hits$Family <- 0
    topCO1hits$Genus <- 0
    topCO1hits$Species <- 0
    topCO1hits$read_hits <- 0
    topCO1hits$percentage_of_CO1_reads <- 0
  }
  
  
}

if ((is.null(CO1microbiome))) { 
  
  topCO1hits <- as.data.frame(matrix(nrow=1, ncol=5))
  colnames(topCO1hits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_CO1_reads")
  
  topCO1hits$Family <- 0
  topCO1hits$Genus <- 0
  topCO1hits$Species <- 0
  topCO1hits$read_hits <- 0
  topCO1hits$percentage_of_CO1_reads <- 0
}




#LSU





LSUmicrobiome <- pavian::read_report(xargs$LSUmicrobiome)


if (!(is.null(LSUmicrobiome))) {
  
  topspecies <- subset(LSUmicrobiome, LSUmicrobiome$taxRank=="S")
  
  topspeciesordered <- topspecies[order(-topspecies$taxonReads),]
  
  if (nrow(topspeciesordered) >=1) {
    
    if (nrow(topspeciesordered) > 10) {
      topspeciestop10 <- topspeciesordered[1:10,]
    }
    if (nrow(topspeciesordered) <= 10) {
      topspeciestop10 <- topspeciesordered
    }
    
    topLSUhits <- as.data.frame(matrix(nrow=nrow(topspeciestop10), ncol=5))
    colnames(topLSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_LSU_reads")
    
    strings <- str_extract(topspeciestop10$taxLineage, "f_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="f_",replacement="")
    topLSUhits$Family <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "g_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="g_",replacement="")
    topLSUhits$Genus <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "s_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="s_",replacement="")
    topLSUhits$Species <- unlist(strings)
    
    topLSUhits$read_hits <- topspeciestop10$taxonreads
    topLSUhits$percentage_of_LSU_reads <- topspeciestop10$percentage
    
  }
  
  if (nrow(topspecies) ==0) { 
    
    topLSUhits <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(topLSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_LSU_reads")
    
    topLSUhits$Family <- 0
    topLSUhits$Genus <- 0
    topLSUhits$Species <- 0
    topLSUhits$read_hits <- 0
    topLSUhits$percentage_of_LSU_reads <- 0
  }
  
  
}

if ((is.null(LSUmicrobiome))) { 
  
  topLSUhits<- as.data.frame(matrix(nrow=1, ncol=5))
  colnames(topLSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_LSU_reads")
  
  topLSUhits$Family <- 0
  topLSUhits$Genus <- 0
  topLSUhits$Species <- 0
  topLSUhits$read_hits <- 0
  topLSUhits$percentage_of_LSU_reads <- 0
}




#SSU
SSUmicrobiome <- pavian::read_report(xargs$SSUmicrobiome)



SSUmicrobiome <- pavian::read_report(xargs$SSUmicrobiome)


if (!(is.null(SSUmicrobiome))) {
  
  topspecies <- subset(SSUmicrobiome, LSUmicrobiome$taxRank=="S")
  
  topspeciesordered <- topspecies[order(-topspecies$taxonReads),]
  
  if (nrow(topspeciesordered) >=1) {
    
    if (nrow(topspeciesordered) > 10) {
      topspeciestop10 <- topspeciesordered[1:10,]
    }
    if (nrow(topspeciesordered) <= 10) {
      topspeciestop10 <- topspeciesordered
    }
    
    topSSUhits <- as.data.frame(matrix(nrow=nrow(topspeciestop10), ncol=5))
    colnames(topSSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_SSU_reads")
    
    strings <- str_extract(topspeciestop10$taxLineage, "f_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="f_",replacement="")
    topSSUhits$Family <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "g_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="g_",replacement="")
    topSSUhits$Genus <- unlist(strings)
    
    strings <- str_extract(topspeciestop10$taxLineage, "s_.*")
    strings <- lapply(strings, gsub,pattern="\\|.*",replacement="")
    strings <- lapply(strings, gsub,pattern="s_",replacement="")
    topSSUhits$Species <- unlist(strings)
    
    topSSUhits$read_hits <- topspeciestop10$taxonreads
    topSSUhits$percentage_of_SSU_reads <- topspeciestop10$percentage
    
  }
  
  if (nrow(topspecies) ==0) { 
    
    topSSUhits <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(topSSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_SSU_reads")
    
    topSSUhits$Family <- 0
    topSSUhits$Genus <- 0
    topSSUhits$Species <- 0
    topSSUhits$read_hits <- 0
    topSSUhits$percentage_of_LSU_reads <- 0
  }
  
  
}

if ((is.null(SSUmicrobiome))) { 
  
  topSSUhits<- as.data.frame(matrix(nrow=1, ncol=5))
  colnames(topSSUhits) <- c("Family", "Genus", "Species", "read_hits", "percentage_of_SSU_reads")
  
  topSSUhits$Family <- 0
  topSSUhits$Genus <- 0
  topSSUhits$Species <- 0
  topSSUhits$read_hits <- 0
  topSSUhits$percentage_of_LSU_reads <- 0
}


topCO1hits$Sample <- NAMES
topLSUhits$Sample <- NAMES
topSSUhits$Sample <- NAMES


write.table(topCO1hits,file=(paste0(outtablespath,NAMES,"_top10returnedspecies_CO1.txt")),sep="\t",row.names=FALSE)
write.table(topLSUhits,file=(paste0(outtablespath,NAMES,"_top10returnedspecies_LSU.txt")),sep="\t",row.names=FALSE)
write.table(topSSUhits,file=(paste0(outtablespath,NAMES,"_top10returnedspecies_SSU.txt")),sep="\t",row.names=FALSE)



# Specialised pathways: Viral identification

# 1. Viral signatures detected, in total, and in unassigned contigs. 

# 2. If ORF annotation undertaken... (?) 


# Combined barplots for assignments of reads start to finish 






# First. Remove extra NAs from diamondraws98%

# Make a 50% for Diamond because I'm calling at superkingdom level
# Can go back and remove or adjust as necessary. Not currently sure what the optimal would be here.

if (dodiamondraws == 'yes') {
  
  Diamondrawshitsident50 <- subset(Diamondrawshits,Diamondrawshits$pident >=50)
  
  if (dohostdetect =="yes") {
    
    if( hostspecies[1] !="NA" ) { 
      
      rows <- grep(hostspecies,Diamondrawshitsident50$species)
      
      if (length(rows)>=1) {
        
        Diamondrawshitsident50nohost <- Diamondrawshitsident50[-rows,]
      }
      if (length(rows)<1) {
        
        Diamondrawshitsident50nohost <- Diamondrawshitsident50
      }
    }
  }
  
  if (dohostdetect =="no") {
    
    Diamondrawshitsident50nohost <- Diamondrawshitsident50
    
  }
  
  idvec <- unique(Diamondrawshitsident50nohost$superkingdom)
  
  
  contigs <- as.data.frame(matrix(nrow=length(idvec), ncol=12))
  colnames(contigs) <- c("Frequency","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident")
  
  
  
  for (i in c(1:length(idvec))) {
    
    contigssubset <- subset(Diamondrawshitsident50nohost,Diamondrawshitsident50nohost$superkingdom == idvec[i])
    
    contigs[i,2:9] <- contigssubset[1,11:18]
    contigs[i,1] <- nrow(contigssubset)
    
  }
  
  freqtablerawdiamond <- contigs
  
  if (dohostdetect =="no") {
    
    allassignedfreqsnohost <- allassignedfreqs
  }
  
  # Repeat for contigs results
  idvec <- unique(allassignedfreqsnohost$superkingdom)
  
  contigs <- as.data.frame(matrix(nrow=length(idvec), ncol=12))
  colnames(contigs) <- c("Frequency","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident")
  
  for (i in c(1:length(idvec))) {
    
    contigssubset<- subset(allassignedfreqsnohost,allassignedfreqsnohost$superkingdom == idvec[i])
    
    contigs[i,2:9] <- contigssubset[1,4:11]
    contigs[i,1] <- sum(contigssubset$freq)
    
  }
  
  freqtablecontigs <- contigs
  
  # merge and final repeat
  allassigned <- rbind(freqtablecontigs,freqtablerawdiamond)
  
  # Need to re add up the families that are multiples again
  
  idvec <- unique(allassigned$superkingdom)
  
  contigs <- as.data.frame(matrix(nrow=length(idvec), ncol=12))
  colnames(contigs) <- c("Frequency","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident")
  
  for (i in c(1:length(idvec))) {
    
    contigssubset<- subset(allassigned,allassigned$superkingdom == idvec[i])
    
    contigs[i,2:12] <- contigssubset[1,2:12]
    contigs[i,1] <- sum(contigssubset$Frequency)
    
  }
  
  
  allassignedordered <- contigs[order(-contigs$Frequency),]
  
  allassignedorderedshort <- as.data.frame(matrix(nrow=nrow(allassignedordered), ncol = 2))
  
  colnames(allassignedorderedshort) <- c("superkingdom","frequency")
  
  allassignedorderedshort$superkingdom <- allassignedordered$superkingdom
  allassignedorderedshort$frequency <- allassignedordered $Frequency
  
  allassignedorderedshort<- allassignedorderedshort[1:4,]
  
  Eukaryotesraws <- subset(allassignedorderedshort,allassignedorderedshort$superkingdom=="Eukaryota")
  
  Bacteriaraws <- subset(allassignedorderedshort,allassignedorderedshort$superkingdom=="Bacteria")
  
  Virusesraws <- subset(allassignedorderedshort,allassignedorderedshort$superkingdom=="Viruses")
  
  Eukaryotesraws[1,2] <- nrow(subset(Diamondrawshitsident50nohost,Diamondrawshitsident50nohost$superkingdom=="Eukaryota"))
  
  Bacteriaraws[1,2] <- nrow(subset(Diamondrawshitsident50nohost,Diamondrawshitsident50nohost$superkingdom=="Bacteria"))
  
  Virusesraws[1,2] <- nrow(subset(Diamondrawshitsident50nohost,Diamondrawshitsident50nohost$superkingdom=="Viruses"))
  
  
}


if (dodiamondraws == 'no') {
  Eukaryotesraws <- as.data.frame(matrix(nrow=1,ncol=2))
  Eukaryotesraws[1,2] <- 0
  Bacteriaraws <- as.data.frame(matrix(nrow=1,ncol=2))
  Bacteriaraws[1,2] <- 0
  Virusesraws <- as.data.frame(matrix(nrow=1,ncol=2))
  Virusesraws[1,2] <- 0
  
}



 


# Generate approx superkingdom read hits across all analyses.

# reads assigned to contigs
Bacteriaspeciessum <- sum(Bacteriaspecies$Frequency)
Eukaryotesspeciessum <- sum(Eukaryotesspecies$Frequency)
Viralspeciessum <- sum(Viralspecies$Frequency)


Eukaryotes <- Eukaryotesraws[1,2] + Eukaryotesspeciessum
Bacteria <- Bacteriaraws[1,2] + Bacteriaspeciessum
Viruses <- Virusesraws[1,2] + Viralspeciessum



# Generate other summ stats results
  

# Values for final bar graph
QC <- (summary_reads_table$filtered_reads)
PhiX <- summary_reads_table$Phix_filtered
CO1 <- summary_reads_table$CO1_filtered
LSU <- summary_reads_table$LSU_filtered
SSU <- summary_reads_table$SSU_filtered
hosthits <- (summary_contigs_table$Raw_reads_assigned_to_host_genome + summary_contigs_table$Raw_reads_assigned_to_host_aligned_contigs + summary_contigs_table$Raw_reads_assigned_to_host_blastn_diamondx)
hosthits <- round(hosthits)
unassignedcontigs <- readsfromunassignedcontigs
unassignedreads <- summary_contigs_table$raw_reads_not_assigned

# create final table for barplot

barplotdf <- as.data.frame(matrix(nrow=1,ncol=11))

colnames(barplotdf) <- c("Low quality reads","PhiX contamination","CO1 hits","LSU hits","SSU hits","Host species","Eukaryotes","Bacteria","Viruses","Reads from unassigned contigs","Unassigned reads")

barplotdf[1,1] <- QC 
barplotdf[1,2] <- PhiX 
barplotdf[1,3] <- CO1 
barplotdf[1,4] <- LSU
barplotdf[1,5] <- SSU 
barplotdf[1,6] <- hosthits
barplotdf[1,7] <- Eukaryotes
barplotdf[1,8] <- Bacteria
barplotdf[1,9] <- Viruses
barplotdf[1,10] <- unassignedcontigs
barplotdf[1,11] <- unassignedreads
barplotdf$Sample <- NAMES

#The reads that partial assign to the host genomes are getting double counted with the raw reads. I need to subtract them off again 


extrareads <- (summary_reads_table$remaining_reads - (hosthits + Eukaryotes + Bacteria + Viruses + unassignedcontigs + unassignedreads ))

unassignedreads2 <- (unassignedreads + extrareads)

barplotdf[1,11] <- unassignedreads2



write.table(barplotdf,file=(paste0(outtablespath,NAMES,"Summary_assignment_reads_for_plot_generation.txt")),sep="\t",row.names=FALSE)

rm(Diamondrawshitshighacc)

save.image(file = paste0(outtablespath,NAMES,"_gather_summary_files_R_environment.Rdata"))


