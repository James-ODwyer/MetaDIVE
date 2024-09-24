# Script for subsetting and extracting valuble information from Diamond blast

# Diamond blast returns the top X (usually set to 10 to 20) hits for each contig. 
# I want to use this information to assign contigs with a putative assignment to species
# or if not possible, to the most recent ancestral grouping
# A quick method for this is to just go by what the top hit is and compile. This has some issues though
# If i took the top hit and the top hit has a vague tax id (e.g., to genus or family predicted) 
# Then that is little information, it might also be that the e values between them are practically
# identical (e.g.,  what is the functional difference betweek 1.5*10^-25 and 1.3*10^25)
# I need a system which can 
# 1. extract the species when there is a clear top hit
# 2. read species Ids across all hits per contig for when all hits are similar
# 3. compile a most useful finding based on what all the contigs say together when the hit rate is 
# practically identical.

# The first question is whether bitscore or evalue provides the optimal measure to use.


#Libraries are scattered throughout as some have common shared function names (e.g., count)


library(plyr)
library(dplyr)
library(argparse)
library("taxonomizr")
library("stringr")
library("foreach")
library("doParallel")
#library("tidyverse")
library(seqinr)
library(tidyr)

#args <- commandArgs(TRUE)
parser <- ArgumentParser(description= 'Informing Diamond blast')

parser$add_argument('--inputdiamond', '-i', help= 'I am the input diamond file')
parser$add_argument('--output', '-o', help= 'Output file for assigned contigs')
parser$add_argument('--programdir', '-p', help= 'working program directory')
parser$add_argument('--savdir', '-s', help= 'working directory for saving contig matches')
#parser$add_argument('--savcontig', '-S', help= 'filename for saving contig matches')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--threads', '-t', help= 'Number of threads')
parser$add_argument('--Log', '-l', help= 'Log of data')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')
parser$add_argument('--clustfile', '-C', help= 'Cluster_file_for_reduced_blastn_reads')



xargs<- parser$parse_args()

# Function to parse in Diamond with built in redundancy (sometimes Diamond/blastn outputs incorrect column numbers because tax info is missing
Diamond_lines<- readLines(xargs$inputdiamond)

# Function to replace '\n' with an actual newline character

replace_newline_within_quotes <- function(line) {
  parts <- strsplit(line, "\"")[[1]]
  for (i in seq_along(parts)) {
    if (i %% 2 == 0) {
      parts[[i]] <- gsub("\n", " ", parts[[i]])
    }
  }
  paste(parts, collapse = "\"")
}

# Function to parse in cluster file for what individual reads compose each clustered read.
parse_clstr <- function(clstr_file) {
  clusters <- list()
  current_cluster <- NULL
  
  lines <- readLines(clstr_file)
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_cluster <- sub("^>Cluster ", "", line)
      clusters[[current_cluster]] <- c()
    } else {
      read_id <- sub(".*>([^.]+)\\.\\.\\..*", "\\1", line)
      clusters[[current_cluster]] <- c(clusters[[current_cluster]], read_id)
    }
  }
  
  clusters_df <- do.call(rbind, lapply(names(clusters), function(cluster) {
    data.frame(cluster_rep = cluster, read_id = clusters[[cluster]], stringsAsFactors = FALSE)
  }))
  
  return(clusters_df)
}

# Parse the .clstr file


# read in cluster file location
clstr_file <- xargs$clustfile
# read in and parse the cluster file
clusters_df <- parse_clstr(clstr_file)




# Apply the function to all lines
modified_lines <- lapply(Diamond_lines, replace_newline_within_quotes)
modified_lines  <- unlist(modified_lines)
# Read the modified lines using read.table
Diamond_output <- read.table(text = modified_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

# Replace the unique character sequence back to '\n' within the dataframe
Diamond_output[] <- lapply(Diamond_output, function(col) gsub("###NEWLINE###", "\n", col))

rm(modified_lines,Diamond_lines)


outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$savdir
outtablespath <- paste0(outtablespath1,outtablespath2)
AccessionNamenode <- xargs$Accnode

NAMES <- xargs$name
n.cores <- xargs$threads


# This is no longer the contigs, this is the raw reads tags for viruses
assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
assigned_contigs <- paste0(assigned_contigs1,assigned_contigs2)

contigsvdir <-xargs$savdir

log <- xargs$Log

# Set to correct directory.
# sink all cats and prints to the snakemake log files 
#sink(log)

cat(paste0("printing results tables to", outtablespath))
cat(paste0("printing matched contigs to", assigned_contigs))

# give proper column names
colnames(Diamond_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")

cat(paste0("Succesfully read in Diamond file for sample", NAMES,"\n"))

cat("Head 10 \n ")
print(Diamond_output[1:10,])

#Extract distinct contigs which generated hits

b <- Sys.time()


# Total time <5min on personal computer therefore not benchmarked


# This section checks to see whether the taxid of the optimal contig hits is one or multiple taxids
# This is important for downstream knowning closest species 
# but also important for actually generating larger taxonomic structure in the next step


Diamond_output$bitscore <- as.numeric(Diamond_output$bitscore)
Diamond_output$pident<- as.numeric(Diamond_output$pident)
Diamond_output$length<- as.numeric(Diamond_output$length)


c <- Sys.time()

# about 10% of returned alignments from Diamond don't have a tax id looking closer these are usually a gene fragment which for whatever reason
# diamond may not be coded to interpret their meta data correctly. 
# To address this as fast as possible I have extracted the species name from stitle and ran that through the unique subsetting and grep sorting below
# to minimise the number of instances where taxonomizr has to search the data base
cat(paste0(" Identification of taxonomy from tax ids  ", "\n"))
cat(paste0(Sys.time(), "\n"))
Diamond_output$staxidreduced<- gsub(pattern = ";.*",replacement = "",x = Diamond_output$staxids)
Diamond_output$sp <- gsub(pattern = ".*\\[",replacement = "",x = Diamond_output$stitle)
Diamond_output$sp <- gsub(pattern = "\\]",replacement = "",x = Diamond_output$sp)

cat(paste0(" Finished extracting species names ", "\n"))

missingidx <- grep("^$",Diamond_output$staxidreduced)

cat(paste0(" Finished generating the missing tax id list ", "\n"))

Diamond_outputmissingonly <- Diamond_output[missingidx,]



uniquespmissing<- dplyr::distinct(Diamond_outputmissingonly, sp, .keep_all = TRUE)

cat(paste0(" Finished subsetting by unique species ", "\n"))

if (nrow(uniquespmissing) >=1) {

taxids <- as.data.frame(matrix(nrow=nrow(Diamond_outputmissingonly),ncol=2))
taxidsunique <- as.data.frame(matrix(nrow=nrow(uniquespmissing),ncol=2))
taxids[,1] <-Diamond_outputmissingonly$sp
taxidsunique[,1] <-uniquespmissing$sp

for ( i in c(1:nrow(taxidsunique))) {
  
  taxidsunique[i,2] <- taxonomizr::getId(uniquespmissing$sp[i], sqlFile=AccessionNamenode)
  
}


cat(paste0(" Finished getID taxonomizr ", "\n"))


cat(paste0(" Dimensions for taxids  ", "\n"))
dim(taxids)

cat(paste0(" Dimensions for taxidsunique ", "\n"))
dim(taxidsunique) 


# slight error in the class assignment in the second column of taxidsunique.
taxidsunique$V2 <- as.character(taxidsunique$V2)

# The grep failed because of regex strings within the text, specifically the . and () 
# I would change the grep to fixed but I need to start/finish characters to minimise 
# returning multiple indexes when for example there are multiple specific strains e.g., E coli
# S1(2020) and just E coli.
# Solution is to remove the characters with gsub before the large grep and recombine

taxids$V1 <- gsub(pattern = ".",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ".",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = "(",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = "(",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = ")",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ")",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = "+",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = "+",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = ":",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ":",replacement = "",x = taxidsunique$V1,fixed = TRUE)

for (i in c(1:nrow(taxids))) {
  
  grep(paste0("^",taxids[i,1],"$"),taxidsunique$V1) -> idxval
  taxids[i,2] <-taxidsunique[idxval,2]
}
cat(paste0(" Finished grep converting unique species back to full taxid list ", "\n"))
cat(paste0(Sys.time(), "\n"))

cat(paste0(" Dimensions for diamond output missing only ", "\n"))
dim(Diamond_outputmissingonly) 

cat(paste0(" number of values in missingidx ", "\n"))
length(missingidx) 

cat(paste0(" Dimensions for diamond output ", "\n"))
dim(Diamond_output) 


Diamond_outputmissingonly$staxidreduced <- taxids[,2]

Diamond_output[missingidx,10] <- Diamond_outputmissingonly$staxidreduced


}


Diamond_output$staxidreduced <- as.numeric(Diamond_output$staxidreduced)

Diamond_output <- Diamond_output[order(Diamond_output$staxidreduced),]

Diamond_output$staxidreduced <- as.character(Diamond_output$staxidreduced)

row.names(Diamond_output) <- 1:nrow(Diamond_output)

contigsassignedunique <- Diamond_output[!duplicated(Diamond_output$staxidreduced), ]

contigsassignedindexes <- as.data.frame(matrix(nrow=nrow(contigsassignedunique),ncol=3))

contigsassignedindexes$V1 <- contigsassignedunique$staxidreduced
contigsassignedindexes$V2 <- row.names(contigsassignedunique)


for ( i in c(1:nrow(contigsassignedindexes))) {
  
  if ( i < nrow(contigsassignedindexes)) {
    
    contigsassignedindexes$V3[i] <- ( as.numeric(contigsassignedindexes$V2[i+1]) - 1 )
    
  }
  
  if ( i == nrow(contigsassignedindexes)) {
    
    contigsassignedindexes$V3[i] <- (as.numeric(nrow(Diamond_output)))
    
  }
  
}

colnames(contigsassignedindexes) <- c("staxidreduced", "startingrow","endingrow")

#
#
#
cat(paste0(" Completed identifying the taxid of unassigned samples using their species names  ", "\n"))

cat(paste0(" Completed identifying the taxid of unassigned samples using their species names  ", "\n"))
cat(paste0(" Starting taxonomy identification from all known taxids ", "\n"))

taxids <- as.data.frame(matrix(nrow=nrow(Diamond_output),ncol=9))
taxidsunique <- as.data.frame(matrix(nrow=nrow(contigsassignedunique),ncol=9))

taxids$V1 <- Diamond_output$staxidreduced
taxidsunique$V1 <- contigsassignedunique$staxidreduced
taxids$V2 <- Diamond_output$sseqid
taxidsunique$V2 <- contigsassignedunique$sseqid

cat(paste0(" Starting taxonomizr getTaxonomy known taxids ", "\n"))


if (nrow(taxidsunique)>0) {
for (i in c(1:nrow(contigsassignedunique))) {
  
  taxidsunique[i,2:8] <- taxonomizr::getTaxonomy(contigsassignedunique$staxidreduced[i], sqlFile=AccessionNamenode)
  
  values <- taxonomizr::getRawTaxonomy(contigsassignedunique$staxidreduced[i], sqlFile=AccessionNamenode)
  
  if (!is.null(values[[1]][1])) {
    values[[1]][1] -> taxidsunique[i,9]
  }
  if (is.null(values[[1]][1])) {
    taxidsunique[i,8] -> taxidsunique[i,9]
  }
      
  }
  
  if (is.na(taxidsunique[i,8])) {
    
    spname <- str_extract_all(contigsassignedunique$stitle[i], "\\[(.*?)\\]")[[1]]
    spname2 <- str_replace_all(spname, "\\[|\\]", "")
    value <- taxonomizr::getId(spname2, sqlFile=AccessionNamenode)
    if (length(value) > 0) {
    taxidsunique[i,2:8] <- taxonomizr::getTaxonomy(value, sqlFile=AccessionNamenode)
    
    
    values <- taxonomizr::getRawTaxonomy(value, sqlFile=AccessionNamenode)
    
  if (!is.null(values[[1]][1])) {
    values[[1]][1] -> taxidsunique[i,9]
  }
  if (is.null(values[[1]][1])) {
    taxidsunique[i,8] -> taxidsunique[i,9]
  }
    
    
  }
  }
  
  
}


d <- Sys.time()

cat(paste0(NAMES," Finished taxid conversion to taxonomy ", "\n"))
cat(paste0(Sys.time(), "\n"))

taxidsunique<- as.data.frame(taxidsunique, ncol=9)
colnames(taxidsunique) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies") 


taxidsunique$startrow <- as.numeric(contigsassignedindexes$startingrow)
taxidsunique$endrow <- as.numeric(contigsassignedindexes$endingrow)


e <- Sys.time()



for ( i in c(1:nrow(taxidsunique))) {
  
  
  taxidsunique$startrow[i] ->start
  taxidsunique$endrow[i] ->fin
  
  
  taxids[start:fin,2:9] <- taxidsunique[i,2:9]
  
}

f <- Sys.time()

cat(paste0(" Finished grep converting unique taxids back to diamond output ", "\n"))
cat(paste0(Sys.time(), "\n"))

colnames(taxids) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies") 

Diamond_output$sp <- NULL

Diamond_output$superkingdom <- taxids$superkingdom
Diamond_output$phylum<- taxids$phylum
Diamond_output$class<- taxids$class
Diamond_output$order<- taxids$order
Diamond_output$family<- taxids$family
Diamond_output$genus<- taxids$genus
Diamond_output$species<- taxids$species
Diamond_output$subspecies<- taxids$subspecies





rm(taxids,taxidsunique,uniquespmissing)

contigsassigned <- Diamond_output



contigsassigned$qseqid <- gsub("\t", " ", contigsassigned$qseqid)
contigsassigned$sseqid <- gsub("\t", " ", contigsassigned$sseqid)
contigsassigned$species <- gsub("\t", " ", contigsassigned$species)
contigsassigned$staxids <- gsub("\t", " ", contigsassigned$staxids)
contigsassigned$staxids <- gsub("/", "-", contigsassigned$staxids)
contigsassigned$stitle <- gsub("[[:punct:]]", "", contigsassigned$stitle)
contigsassigned$species <- gsub("[[:punct:]]", "", contigsassigned$species)
contigsassigned$subspecies <- gsub("\t", " ", contigsassigned$subspecies)
contigsassigned$subspecies <- gsub("[[:punct:]]", "", contigsassigned$subspecies)

# map back cluster file and expand the contigsassignment blastn results

#This was really slow for standard for loop. The following speed ups are implemented
# Rbind occurrs once after the loop, intermitten results are stored in a list.
# dplyr filters are used instead of subset
# After each iteration of i, all assigned values are removed from the cluster_df table to speed up subsequent iterations and stop it from searching
# the full file every time. 



results_list <- vector("list", nrow(contigsassigned))

# Initialize a counter for timing
iteration_counter <- 0

# Create a copy of clusters_df for searching
clusters_df_search <- clusters_df

# Main for loop
for (i in seq_len(nrow(contigsassigned))) {
  
  # Measure time for every 500 iterations
  if (i %% 500 == 1) {
    iteration_start_time <- Sys.time()
  }
  
  # Perform the search on the copy
  index <- grep(contigsassigned$qseqid[i], clusters_df_search$read_id)
  
  if (length(index) > 0) {
    index2 <- as.numeric(clusters_df_search$cluster_rep[index])
    
    clusters_dfsubset <- clusters_df_search %>%
      filter(cluster_rep %in% index2)
    
    rows_to_add <- nrow(clusters_dfsubset)
    
    blastn_readsubset <- data.frame(matrix(nrow = rows_to_add, ncol = ncol(contigsassigned)))
    colnames(blastn_readsubset) <- colnames(contigsassigned)
    
    blastn_readsubset[, 1] <- clusters_dfsubset$read_id
    blastn_readsubset[, 2:18] <- contigsassigned[i, 2:18]
    
    # Store the result in the list
    results_list[[i]] <- blastn_readsubset
    
    # Remove the identified rows from the search copy
    clusters_df_search <- clusters_df_search %>%
      filter(!cluster_rep %in% index2)
  }
  
  # Increment the counter
  iteration_counter <- iteration_counter + 1
  
  # Print the time taken for the last 500 iterations
  if (i %% 500 == 0) {
    iteration_end_time <- Sys.time()
    iteration_time_taken <- iteration_end_time - iteration_start_time
    cat("Time taken for iterations", (i-499), "to", i, ":", iteration_time_taken, "\n")
  }
}

# Combine all results once after the loop
contigsassigned_extended <- do.call(rbind, results_list)

# Remove NULL elements from the list (if any)
contigsassigned_extended <- contigsassigned_extended[!sapply(contigsassigned_extended, is.null), ]

# Convert the list to a data frame
contigsassigned_extended <- bind_rows(contigsassigned_extended)














contigsassigned_extended <- contigsassigned_extended %>%
  distinct()

rm(Diamond_output)

cat(paste0(NAMES," Finished reindexing of taxa results ", "\n"))
cat(paste0(Sys.time(), "\n"))

contigsassigned_extended$stitle = substr(contigsassigned_extended$stitle,1,50)
# generate summary stats for blast

# Now this section is used to update the original assembled contig names so that they include
# the additional species information (This step may be superfluous later on but it may also be useful
# to compare to see how genome binning works compared to direct assignment)

g <- Sys.time()

i=1
j=1

# Revised index and order based code here
# It works first by extracting the contig number from the contigs assignedworking dataframe (the data frame containing blast hits plus taxonomy data)
# Then the df is ordered by the contig number
# Then extract the names from the megacontigs file (from Megahit) 
# Then extract the contig number from the names as a new column [,2]
# then creating a basic index for what row each contig is from [,3]
# Then filtering the contigs in the Megahit names to only include the rows which have contigs which were blast assigned. 
# Then copy over all name information from the full megahit file into a final names results  (it needs to be the full one to not lead to Na errors)
#then it is a simple loop for each row in the megahit contigs hits so that the full data from contigs assigned is put into the name of each correct row
# (row based on the row index created in the Megahit contigs) and finally reading that full final names which has some regular names and some overridden names
# at the correct spot straight into the contigs information and saving the fasta. 

cat(paste0(NAMES," Starting summary statistics generation per superkingdom ", "\n"))
cat(paste0(Sys.time(), "\n"))

Viruses <- subset(contigsassigned_extended, contigsassigned_extended$superkingdom=="Viruses")


totalcontigsassigned <- nrow(contigsassigned_extended)


Virusessinglesp <- Viruses
Virusessinglespspcounts <-plyr::count(Virusessinglesp$species)
nrow(Virusessinglesp)
sum(Virusessinglesp$species!="NA")


if (nrow(Virusessinglespspcounts) >0) {
  
  Virusessinglespspcounts$V3 <- NA 
  
  for (i in c(1:nrow(Virusessinglespspcounts))) {
    
    Virusessinglespspcounts$V3[i] <- ((Virusessinglespspcounts$freq[i]/totalcontigsassigned)*100)
    
  }
  
  
  colnames(Virusessinglespspcounts) <- c("Species", "Frequency", "Percentage_total_assigned")
  Virusessinglespspcounts<- Virusessinglespspcounts[order(Virusessinglespspcounts$Frequency, decreasing=TRUE),]
  
}


superkingdom <- plyr::count(contigsassigned_extended$superkingdom)

cat(paste0(NAMES," Finished summary statistics generation per superkingdom ", "\n"))
cat(paste0(Sys.time(), "\n"))





viralreadsnames <- Viruses$qseqid




write.table(Virusessinglespspcounts, file=(paste0(outtablespath,NAMES,"_Virus_species_hits_summary.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Virusessinglespspcounts),sep="\t")


# Write the table for the raw read names of viruses

write.table(viralreadsnames ,file=assigned_contigs,sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)

write.table(contigsassigned_extended, file=(paste0(outtablespath,NAMES,"_all_reads_assignments.txt")),sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)



cat(paste0("Finished individual ",NAMES,"\n"))
cat(paste0(Sys.time(), "\n"))




