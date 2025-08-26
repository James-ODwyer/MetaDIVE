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

Diamond_output$bitscore <- as.numeric(Diamond_output$bitscore)
Diamond_output$length   <- as.numeric(Diamond_output$length)
Diamond_output$pident   <- as.numeric(Diamond_output$pident)

# --- Coalesce multiple hits per qseqid into: primary + top-3 alternates (genus/species), same code as original diamond contigs
#take first taxid if semicolon-delimited
reduce_taxid <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  sub(";.*", "", x, perl = TRUE)
}

# Resolve genus/species for a taxid (for alternates only)
get_genus_species_from_taxid <- function(taxid, sql) {
  if (is.na(taxid) || taxid == "") return(c(genus = "NONE", species = "NONE"))
  suppressWarnings({
    tx <- try(taxonomizr::getTaxonomy(taxid, sqlFile = sql), silent = TRUE)
  })
  if (inherits(tx, "try-error") || is.null(tx)) return(c(genus = "NONE", species = "NONE"))
  g <- if (!is.null(tx[,"genus"])   && !is.na(tx[,"genus"]))   tx[,"genus"]   else "NONE"
  s <- if (!is.null(tx[,"species"]) && !is.na(tx[,"species"])) tx[,"species"] else "NONE"
  c(genus = as.character(g), species = as.character(s))
}

ALT_WITHIN <- 0.90  # consider alternates within 90% of top bitscore

# Stable order: by read then by best hits
Diamond_output <- Diamond_output[order(Diamond_output$qseqid, -Diamond_output$bitscore, Diamond_output$evalue), ]
qids <- unique(Diamond_output$qseqid)

# Result: one-row-per-read with 6 new columns at the end
coalesced <- as.data.frame(matrix(nrow = length(qids), ncol = ncol(Diamond_output) + 6))
colnames(coalesced) <- c(colnames(Diamond_output),
                         "alternate_genus1","alternate_species1",
                         "alternate_genus2","alternate_species2",
                         "alternate_genus3","alternate_species3")

for (i in seq_along(qids)) {
  wc <- Diamond_output[Diamond_output$qseqid == qids[i], , drop = FALSE]
  if (nrow(wc) == 0) next
  primary <- wc[1, , drop = FALSE]
  
  # candidates within the threshold (incl. primary), dedupe alternates by reduced taxid
  keep_idx <- which(wc$bitscore >= ALT_WITHIN * primary$bitscore[1])
  cand <- wc[keep_idx, , drop = FALSE]
  
  primary_taxid <- reduce_taxid(primary$staxids[1])
  cand$taxid_reduced <- vapply(cand$staxids, reduce_taxid, character(1))
  alt_cand <- cand[!is.na(cand$taxid_reduced) & cand$taxid_reduced != primary_taxid, , drop = FALSE]
  alt_cand <- alt_cand[!duplicated(alt_cand$taxid_reduced), , drop = FALSE]
  alt_taxids <- head(alt_cand$taxid_reduced, 3)
  
  # resolve genus/species for up to 3 alternates (cheap: 0–3 calls per read)
  alt_gs <- lapply(alt_taxids, function(tid) get_genus_species_from_taxid(tid, AccessionNamenode))
  while (length(alt_gs) < 3L) alt_gs <- c(alt_gs, list(c(genus = "NONE", species = "NONE")))
  
  outrow <- primary
  outrow$alternate_genus1    <- alt_gs[[1]]["genus"]
  outrow$alternate_species1  <- alt_gs[[1]]["species"]
  outrow$alternate_genus2    <- alt_gs[[2]]["genus"]
  outrow$alternate_species2  <- alt_gs[[2]]["species"]
  outrow$alternate_genus3    <- alt_gs[[3]]["genus"]
  outrow$alternate_species3  <- alt_gs[[3]]["species"]
  
  coalesced[i, ] <- outrow[1, colnames(coalesced)]
}

# 3) rename new coalesced table back to Diamond_output (downstream expects this name)
Diamond_output <- coalesced
rm(coalesced)

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

# Primary-only table now; still recover missing taxids from stitle
Diamond_output$staxidreduced <- sub(";.*", "", Diamond_output$staxids)
Diamond_output$sp <- sub("\\]", "", sub(".*\\[", "", Diamond_output$stitle))

missingidx <- which(is.na(Diamond_output$staxidreduced) | Diamond_output$staxidreduced == "")
Diamond_outputmissingonly <- Diamond_output[missingidx, , drop = FALSE]

if (nrow(Diamond_outputmissingonly) >= 1) {
  uniquespmissing <- dplyr::distinct(Diamond_outputmissingonly, sp, .keep_all = FALSE)
  
  taxids <- data.frame(Vorig = Diamond_outputmissingonly$sp, V1 = Diamond_outputmissingonly$sp,
                       V2 = NA_character_, stringsAsFactors = FALSE)
  taxidsunique <- data.frame(Vorig = uniquespmissing$sp, V1 = uniquespmissing$sp,
                             V2 = NA_character_, stringsAsFactors = FALSE)
  
  # sanitize ONLY the join keys (V1) for anchored grep
  scrub <- function(x) {
    x <- gsub("\\.", "", x, fixed = TRUE)
    x <- gsub("\\(", "", x, fixed = TRUE)
    x <- gsub("\\)", "", x, fixed = TRUE)
    x <- gsub("\\+", "", x, fixed = TRUE)
    x <- gsub(":",  "", x, fixed = TRUE)
    x
  }
  taxids$V1       <- scrub(taxids$V1)
  taxidsunique$V1 <- scrub(taxidsunique$V1)
  
  # getId on the ORIGINAL (unscrubbed) names
  for (i in seq_len(nrow(taxidsunique))) {
    val <- try(taxonomizr::getId(taxidsunique$Vorig[i], sqlFile = AccessionNamenode), silent = TRUE)
    if (!inherits(val, "try-error") && length(val) > 0) taxidsunique$V2[i] <- as.character(val[1])
  }
  
  # map back scrubbed->id
  for (i in seq_len(nrow(taxids))) {
    idxval <- grep(paste0("^", taxids$V1[i], "$"), taxidsunique$V1)
    if (length(idxval) == 1) taxids$V2[i] <- taxidsunique$V2[idxval]
  }
  Diamond_output$staxidreduced[missingidx] <- taxids$V2
}


#
#
#
cat(paste0(" Completed identifying the taxid of unassigned samples using their species names  ", "\n"))

cat(paste0(" Completed identifying the taxid of unassigned samples using their species names  ", "\n"))
cat(paste0(" Starting taxonomy identification from all known taxids ", "\n"))

# One row per read now → one taxonomy lookup each
taxids <- data.frame(staxidreduced = Diamond_output$staxidreduced,
                     sseqid = Diamond_output$sseqid,
                     stringsAsFactors = FALSE)

taxmap <- matrix(NA_character_, nrow = nrow(Diamond_output), ncol = 8)
colnames(taxmap) <- c("superkingdom","phylum","class","order","family","genus","species","subspecies")

for (i in seq_len(nrow(Diamond_output))) {
  tid <- Diamond_output$staxidreduced[i]
  if (!is.na(tid) && tid != "") {
    tx <- try(taxonomizr::getTaxonomy(tid, sqlFile = AccessionNamenode), silent = TRUE)
    if (!inherits(tx, "try-error") && !is.null(tx)) {
      taxmap[i, ] <- c(ifelse(is.na(tx[,"superkingdom"]), NA, as.character(tx[,"superkingdom"])),
                       ifelse(is.na(tx[,"phylum"]),       NA, as.character(tx[,"phylum"])),
                       ifelse(is.na(tx[,"class"]),        NA, as.character(tx[,"class"])),
                       ifelse(is.na(tx[,"order"]),        NA, as.character(tx[,"order"])),
                       ifelse(is.na(tx[,"family"]),       NA, as.character(tx[,"family"])),
                       ifelse(is.na(tx[,"genus"]),        NA, as.character(tx[,"genus"])),
                       ifelse(is.na(tx[,"species"]),      NA, as.character(tx[,"species"])),
                       ifelse(is.na(tx[,"subspecies"]),   NA, as.character(tx[,"subspecies"])))
    } else {
      # fallback: try species name from stitle
      spname <- sub("\\]", "", sub(".*\\[", "", Diamond_output$stitle[i]))
      val <- try(taxonomizr::getId(spname, sqlFile = AccessionNamenode), silent = TRUE)
      if (!inherits(val, "try-error") && length(val) > 0) {
        tx2 <- try(taxonomizr::getTaxonomy(val[1], sqlFile = AccessionNamenode), silent = TRUE)
        if (!inherits(tx2, "try-error") && !is.null(tx2)) {
          taxmap[i, ] <- c(ifelse(is.na(tx2[,"superkingdom"]), NA, as.character(tx2[,"superkingdom"])),
                           ifelse(is.na(tx2[,"phylum"]),       NA, as.character(tx2[,"phylum"])),
                           ifelse(is.na(tx2[,"class"]),        NA, as.character(tx2[,"class"])),
                           ifelse(is.na(tx2[,"order"]),        NA, as.character(tx2[,"order"])),
                           ifelse(is.na(tx2[,"family"]),       NA, as.character(tx2[,"family"])),
                           ifelse(is.na(tx2[,"genus"]),        NA, as.character(tx2[,"genus"])),
                           ifelse(is.na(tx2[,"species"]),      NA, as.character(tx2[,"species"])),
                           ifelse(is.na(tx2[,"subspecies"]),   NA, as.character(tx2[,"subspecies"])))
        }
      }
    }
  }
}

# project taxonomy back to Diamond_output
Diamond_output$superkingdom <- taxmap[, "superkingdom"]
Diamond_output$phylum       <- taxmap[, "phylum"]
Diamond_output$class        <- taxmap[, "class"]
Diamond_output$order        <- taxmap[, "order"]
Diamond_output$family       <- taxmap[, "family"]
Diamond_output$genus        <- taxmap[, "genus"]
Diamond_output$species      <- taxmap[, "species"]
Diamond_output$subspecies   <- taxmap[, "subspecies"]
rm(taxmap)


Diamond_output$sp <- NULL
# Ensure alternate_* columns sit after all other columns
alt_cols <- c("alternate_genus1","alternate_species1",
              "alternate_genus2","alternate_species2",
              "alternate_genus3","alternate_species3")
ord <- c(setdiff(colnames(Diamond_output), alt_cols), alt_cols)
Diamond_output <- Diamond_output[, ord, drop = FALSE]


#rm(taxids,taxidsunique,uniquespmissing)

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
  index <- grep(paste0(contigsassigned$qseqid[i],"$"), clusters_df_search$read_id)
  
  if (length(index) > 0) {
    index2 <- as.numeric(clusters_df_search$cluster_rep[index])
    
    clusters_dfsubset <- clusters_df_search %>%
      filter(cluster_rep %in% index2)
    
    rows_to_add <- nrow(clusters_dfsubset)
    
    blastn_readsubset <- data.frame(matrix(nrow = rows_to_add, ncol = ncol(contigsassigned)))
    colnames(blastn_readsubset) <- colnames(contigsassigned)
    
    blastn_readsubset[, 1] <- clusters_dfsubset$read_id
blastn_readsubset[, 2:ncol(contigsassigned)] <- contigsassigned[i, 2:ncol(contigsassigned)]
    
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

# Combine all results once after the loop (robust to NULLs now)
results_list <- Filter(Negate(is.null), results_list)
if (length(results_list)) {
  contigsassigned_extended <- bind_rows(results_list)
} else {
  contigsassigned_extended <- contigsassigned[0, ]  # empty, same cols
}













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




