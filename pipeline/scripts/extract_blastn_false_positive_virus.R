# Script for subsetting and extracting valuble information from Blastn

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

# need to decide on a useful pident filter to for these hits. Is a 50% ident to mosquito species valuble? Lax/strict filtering?


library(plyr)
library(dplyr)
library(argparse)
library("taxonomizr")
library("stringr")
library("foreach")
library("doParallel")
#library("tidyverse")
library(seqinr)
library(phylotools)

#args <- commandArgs(TRUE)
parser <- ArgumentParser(description= 'Informing Blastn')

parser$add_argument('--inputblastn', '-i', help= 'I am the input diamond file')
parser$add_argument('--abundances', '-a', help= 'Output file path for abundance hit tables')
parser$add_argument('--programdir ', '-p', help= 'working program directory')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--threads', '-t', help= 'Number of threads')
parser$add_argument('--Log', '-l', help= 'Log of data')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')




xargs<- parser$parse_args()


blastn_lines<- readLines(xargs$inputblastn)


replace_newline_within_quotes <- function(line) {
  parts <- strsplit(line, "\"")[[1]]
  for (i in seq_along(parts)) {
    if (i %% 2 == 0) {
      parts[[i]] <- gsub("\n", " ", parts[[i]])
    }
  }
  paste(parts, collapse = "\"")
}


# Apply the function to all lines
modified_lines <- lapply(blastn_lines, replace_newline_within_quotes)
modified_lines  <- unlist(modified_lines)
modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
# Read the modified lines using read.table
blastn_output <- read.table(text = modified_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

# Replace the unique character sequence back to '\n' within the dataframe
blastn_output[] <- lapply(blastn_output, function(col) gsub("###NEWLINE###", "\n", col))

rm(modified_lines,blastn_lines)



colnames(blastn_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")






outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$abundances
outtablespath <- paste0(outtablespath1,outtablespath2)

AccessionNamenode <- xargs$Accnode

NAMES <- xargs$name
n.cores <- xargs$threads


log <- xargs$Log




#multitaxcluster <- parallel::makeCluster(
#  n.cores, 
#  type = "PSOCK"
#  )
#print(multitaxcluster)

#doParallel::registerDoParallel(cl = multitaxcluster)
#foreach::getDoParRegistered()
#foreach::getDoParWorkers()

# Set to correct directory.
# sink all cats and prints to the snakemake log files 
 sink(log)

if (!file.exists(outtablespath)){

dir.create(outtablespath, recursive = FALSE, mode = "0777")
}

cat(paste0("printing results tables to", outtablespath))



  # give proper column names
  colnames(blastn_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")
  
blastn_output$bitscore <- as.numeric(blastn_output$bitscore)
blastn_output$length <- as.numeric(blastn_output$length)
blastn_output$pident <- as.numeric(blastn_output$pident)

  cat(paste0("Succesfully read in Blastn file for sample", NAMES,"\n"))
  
  cat("Head 10 \n ")
  print(blastn_output[1:10,])
  
  #Extract distinct contigs which generated hits
  unique(blastn_output$qseqid) -> contigs
  
  # General storage matrix (needs to be big enough to store all similar bitvalue hits) I set it 
  # to 20 here as 20 was the number of hits returned
  namesorter <- matrix(nrow= 20,ncol = 1)
  
  #General storage list
  namehits <- list()
  
  # generate data frame to hold finished contigs and assignments
  contigsassigned <- as.data.frame(matrix(nrow = length(contigs), ncol = ncol(Diamond_output) + 6))
  colnames(contigsassigned) <- c(colnames(Diamond_output),
                                 "alternate_genus1","alternate_species1",
                                 "alternate_genus2","alternate_species2",
                                 "alternate_genus3","alternate_species3")
  
  # Helper: reduce semicolon-delimited taxids -> first taxid
  reduce_taxid <- function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    sub(";.*", "", x, perl = TRUE)
  }
  
  get_genus_species_from_taxid <- function(taxid, sql) {
    if (is.na(taxid) || taxid == "") return(c(genus = "NONE", species = "NONE"))
    suppressWarnings({
      tx <- try(taxonomizr::getTaxonomy(taxid, sqlFile = sql), silent = TRUE)
    })
    if (inherits(tx, "try-error") || is.null(tx)) return(c(genus = "NONE", species = "NONE"))
    # Make sure we can index by column name
    g <- if (!is.null(tx[,"genus"])   && !is.na(tx[,"genus"]))   tx[,"genus"]   else "NONE"
    s <- if (!is.null(tx[,"species"]) && !is.na(tx[,"species"])) tx[,"species"] else "NONE"
    c(genus = as.character(g), species = as.character(s))
  }
  
  # Numeric coercions (already done above, safe to repeat)
  Diamond_output$bitscore <- as.numeric(Diamond_output$bitscore)
  Diamond_output$pident   <- as.numeric(Diamond_output$pident)
  Diamond_output$length   <- as.numeric(Diamond_output$length)
  
  ALT_WITHIN <- 0.90  # keep hits with bitscore >= 90% of top hit
  
  a <- Sys.time()
  for (i in seq_along(contigs)) {
    
    wcA <- subset(Diamond_output, Diamond_output$qseqid == contigs[i])
    if (nrow(wcA) == 0) next
    
    # sort by bitscore desc (stable within equal bitscore)
    wc <- wcA[order(-wcA$bitscore), ]
    
    # Primary (top) hit row
    primary <- wc[1, , drop = FALSE]
    
    # Candidates within 10% of top bitscore (including the primary itself)
    keep_idx <- which(wc$bitscore >= ALT_WITHIN * primary$bitscore[1])
    cand <- wc[keep_idx, , drop = FALSE]
    
    # Primary reduced taxid
    primary_taxid <- reduce_taxid(primary$staxids[1])
    
    # Deduplicate candidates by reduced taxid, preserve order by bitscore
    cand$taxid_reduced <- vapply(cand$staxids, reduce_taxid, character(1))
    # Exclude the primary taxid for alternates
    alt_cand <- cand[!is.na(cand$taxid_reduced) & cand$taxid_reduced != primary_taxid, , drop = FALSE]
    alt_cand <- alt_cand[!duplicated(alt_cand$taxid_reduced), , drop = FALSE]
    
    # Take up to 3 alternates (highest bitscore already by order)
    alt_taxids <- head(alt_cand$taxid_reduced, 3)
    
    # Resolve genus/species for alternates
    alt_gs <- lapply(alt_taxids, function(tid) get_genus_species_from_taxid(tid, AccessionNamenode))
    # Pad to length 3 with NONEs if needed
    while (length(alt_gs) < 3L) alt_gs <- c(alt_gs, list(c(genus = "NONE", species = "NONE")))
    
    # Fill output row = primary hit + 6 new columns
    outrow <- primary
    outrow$alternate_genus1  <- alt_gs[[1]]["genus"]
    outrow$alternate_species1<- alt_gs[[1]]["species"]
    outrow$alternate_genus2  <- alt_gs[[2]]["genus"]
    outrow$alternate_species2<- alt_gs[[2]]["species"]
    outrow$alternate_genus3  <- alt_gs[[3]]["genus"]
    outrow$alternate_species3<- alt_gs[[3]]["species"]
    
    # Assign into master (one row per contig)
    contigsassigned[i, ] <- outrow[1, colnames(contigsassigned)]
  }
  b <- Sys.time()
  
  # Total time <5min on personal computer therefore not benchmarked
  
  cat(paste0(NAMES," Completed optimal blast hit identification ", "\n"))
  cat(paste0(Sys.time(), "\n"))
  
  
  # This section checks to see whether the taxid of the optimal contig hits is one or multiple taxids
  # This is important for downstream knowning closest species 
  # but also important for actually generating larger taxonomic structure in the next step
  
  contigsassigned$multiplesp <- NA
  contigsassigned$staxidreduced <- NA
  
  for (i in c(1:nrow(contigsassigned))) {
    
    if (grepl(pattern= ";" , contigsassigned$staxids[i])) {
      contigsassigned$multiplesp[i] <-"yes"
      
      contigsassigned$staxidreduced[i] <- gsub(";.*","", x = contigsassigned$staxids[i], perl = TRUE)
    }
    
    else{
      contigsassigned$multiplesp[i] <-"no"
      
      contigsassigned$staxidreduced[i]<- contigsassigned$staxids[i]
      
    }
    
  }

  
  
  cat(paste0(nrow(contigsassigned),"contigs assigned to reference sequences: ","\n"))
  
c <- Sys.time()

contigsassignedunique<- dplyr::distinct(contigsassigned, staxidreduced, .keep_all = TRUE)

taxids <- as.data.frame(matrix(nrow=nrow(contigsassigned),ncol=9))
taxidsunique <- as.data.frame(matrix(nrow=nrow(contigsassignedunique),ncol=9))
taxids[,1] <-contigsassigned$staxidreduced
taxidsunique[,1] <-contigsassignedunique$staxidreduced

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

d <- Sys.time()

	cat(paste0(NAMES," Finished taxid conversion to taxonomy ", "\n"))
	cat(paste0(Sys.time(), "\n"))
	cat(paste0(" Time taken: ", difftime(d,c), "\n"))




taxidsunique<- as.data.frame(taxidsunique, ncol=9)


e <- Sys.time()

taxids$V1 <- gsub(pattern = ".",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ".",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = "(",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = "(",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = ")",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ")",replacement = "",x = taxidsunique$V1,fixed = TRUE)

colnames(taxidsunique) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies") 

for (i in c(1:nrow(taxids))) {

if (!is.na(taxids[i,1])) {

grep(paste0("^",taxids[i,1],"$"),taxidsunique$staxidreduced) -> idxval

taxids[i,] <-taxidsunique[idxval,]
}

if (is.na(taxids[i,1])) {

taxids[i,2:8] <- NA

}

}

f <- Sys.time()
colnames(taxids) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies") 

contigsassigned$superkingdom <- taxids$superkingdom
contigsassigned$phylum<- taxids$phylum
contigsassigned$class<- taxids$class
contigsassigned$order<- taxids$order
contigsassigned$family<- taxids$family
contigsassigned$genus<- taxids$genus
contigsassigned$species<- taxids$species
contigsassigned$subspecies<- taxids$subspecies


contigsassigned$stitle = substr(contigsassigned$stitle,1,50)


      cat(paste0(NAMES," Finished reindexing of taxa results ", "\n"))
      cat(paste0(Sys.time(), "\n"))

  
# generate summary stats for blast
  
  # Now this section is used to update the original assembled contig names so that they include
  # the additional species information (This step may be superfluous later on but it may also be useful
  # to compare to see how genome binning works compared to direct assignment)
  
      
      # Move alternate_* columns to the end of contigsassigned
      alt_cols <- c("alternate_genus1", "alternate_species1",
                    "alternate_genus2", "alternate_species2",
                    "alternate_genus3", "alternate_species3")
      
      # Ensure these columns exist before reordering
      alt_cols <- alt_cols[alt_cols %in% colnames(contigsassigned)]
      
      
      
g <- Sys.time()
    
  i=1
  j=1
  
    
  contigsassignedworking <- contigsassigned


  # the old script added in a name for contigs already identified. Here I have taken a different approach for the Trinity . i've also deleted the megahit
# and am running them both through as just saving the fasta of those which didnt return a hit here. If i end up needing the positive hits later, a simple
# filter against the megahitdf2 for the base will work. 
  
  contigsassignedworking <- contigsassigned
 



Eukaryotes <- subset(contigsassigned, contigsassigned$superkingdom=="Eukaryota")
Bacteria <- subset(contigsassigned, contigsassigned$superkingdom=="Bacteria")
Viruses <- subset(contigsassigned, contigsassigned$superkingdom=="Viruses")



totalcontigsassigned <- nrow(contigsassigned)

contigsassignedsinglecontig <- contigsassignedunique<- dplyr::distinct(contigsassigned, sseqid, .keep_all = TRUE)

Eukaryotesmultisp <- subset(Eukaryotes, Eukaryotes$multiplesp=="yes")
nrow(Eukaryotesmultisp)
sum(Eukaryotesmultisp$genus!="NA")
Eukaryotesmultispgenuscounts <-plyr::count(Eukaryotesmultisp$genus)
Eukaryotessinglesp <- subset(Eukaryotes, Eukaryotes$multiplesp=="no")
Eukaryotessinglespspcounts <-plyr::count(Eukaryotessinglesp$species)
nrow(Eukaryotessinglesp)
sum(Eukaryotessinglesp$species!="NA")





if (nrow(Eukaryotessinglespspcounts) >0) {

  Eukaryotessinglespspcounts$V3 <- NA  
  
  for (i in c(1:nrow(Eukaryotessinglespspcounts))) {
    
    Eukaryotessinglespspcounts$V3[i] <- ((Eukaryotessinglespspcounts$freq[i]/totalcontigsassigned)*100)
  }
  colnames(Eukaryotessinglespspcounts) <- c("Species", "Frequency", "Percentage_total_assigned")

  Eukaryotessinglespspcounts <- Eukaryotessinglespspcounts[order(Eukaryotessinglespspcounts$Frequency, decreasing=TRUE),]

}
  

if (nrow(Eukaryotesmultispgenuscounts) >0) {


Eukaryotesmultispgenuscounts$V3 <- NA

for (i in c(1:nrow(Eukaryotesmultispgenuscounts))) {
  
  Eukaryotesmultispgenuscounts$V3[i] <- ((Eukaryotesmultispgenuscounts$freq[i]/totalcontigsassigned)*100) 
}

colnames(Eukaryotesmultispgenuscounts ) <- c("Genus", "Frequency", "Percentage_total_assigned")

Eukaryotesmultispgenuscounts <- Eukaryotesmultispgenuscounts[order(Eukaryotesmultispgenuscounts$Frequency, decreasing=TRUE),]


}




Bacteriamultisp <- subset(Bacteria, Bacteria$multiplesp=="yes")
nrow(Bacteriamultisp)
sum(Bacteriamultisp$genus!="NA")
Bacteriamultispgenuscounts <-plyr::count(Bacteriamultisp$genus)
Bacteriasinglesp <- subset(Bacteria, Bacteria$multiplesp=="no")
Bacteriasinglespspcounts <-plyr::count(Bacteriasinglesp$species)
nrow(Bacteriasinglesp)
sum(Bacteriasinglesp$species!="NA")



if (nrow(Bacteriasinglespspcounts) >0) {

  Bacteriasinglespspcounts$V3 <- NA
  
  
  for (i in c(1:nrow(Bacteriasinglespspcounts))) {
    
    Bacteriasinglespspcounts$V3[i] <- ((Bacteriasinglespspcounts$freq[i]/totalcontigsassigned)*100)
  }
  
  
  colnames(Bacteriasinglespspcounts) <- c("Species", "Frequency", "Percentage_total_assigned")
  Bacteriasinglespspcounts<- Bacteriasinglespspcounts[order(Bacteriasinglespspcounts$Frequency, decreasing=TRUE),]
}



if (nrow(Bacteriamultispgenuscounts) >0) {

  
  Bacteriamultispgenuscounts$V3 <- NA
  
for (i in c(1:nrow(Bacteriamultispgenuscounts))) {
  
  Bacteriamultispgenuscounts$V3[i] <- ((Bacteriamultispgenuscounts$freq[i]/totalcontigsassigned)*100) 
}


  
  colnames(Bacteriamultispgenuscounts) <- c("Genus", "Frequency", "Percentage_total_assigned")
  
  Bacteriamultispgenuscounts<- Bacteriamultispgenuscounts[order(Bacteriamultispgenuscounts$Frequency, decreasing=TRUE),]

  }



Virusesamultisp <- subset(Viruses, Viruses$multiplesp=="yes")
nrow(Virusesamultisp)
sum(Virusesamultisp$genus!="NA")
Virusesamultispgenuscounts <-plyr::count(Virusesamultisp$genus)
Virusessinglesp <- subset(Viruses, Viruses$multiplesp=="no")
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


if (nrow(Virusesamultispgenuscounts) >0) {
  
  
  Virusesamultispgenuscounts$V3 <- NA
  
  for (i in c(1:nrow(Virusesamultispgenuscounts))) {
    
    Virusesamultispgenuscounts$V3[i] <- ((Virusesamultispgenuscounts$freq[i]/totalcontigsassigned)*100)
  }
  
  
  colnames(Virusesamultispgenuscounts) <- c("Genus", "Frequency","Percentage_total_assigned")
  Virusesamultispgenuscounts<- Virusesamultispgenuscounts[order(Virusesamultispgenuscounts$Frequency, decreasing=TRUE),]
}




superkingdom <- plyr::count(contigsassigned$superkingdom)



# Generate contig list for host species
# 1. can add the species id generated from CO1-SSU in the script the be imported and use that but for the purposes of full analysis, It is very clear what is the host species
# approx 80-90% of contigs are host! 
# This is after host filtering of raw reads. 

# 





  write.table(superkingdom, file=(paste0(outtablespath,NAMES,"_superkingdoms.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(superkingdom),sep="\t")
  write.table(Eukaryotessinglespspcounts, file=(paste0(outtablespath,NAMES,"_Eukaryotes_sp_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Eukaryotessinglespspcounts),sep="\t")
  write.table(Eukaryotesmultispgenuscounts, file=(paste0(outtablespath,NAMES,"_Eukaryotes_genus_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Eukaryotesmultispgenuscounts),sep="\t")
  write.table(Bacteriasinglespspcounts, file=(paste0(outtablespath,NAMES,"_Bacteria_species_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Bacteriasinglespspcounts),sep="\t")
  write.table(Bacteriamultispgenuscounts, file=(paste0(outtablespath,NAMES,"_Bacteria_genus_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Bacteriamultispgenuscounts),sep="\t")
  write.table(Virusessinglespspcounts, file=(paste0(outtablespath,NAMES,"_Virus_species_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Virusessinglespspcounts),sep="\t")
  write.table(Virusesamultispgenuscounts, file=(paste0(outtablespath,NAMES,"_Virus_genus_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Virusesamultispgenuscounts),sep="\t")

  write.table(contigsassigned, file=(paste0(outtablespath,NAMES,"_Contigsallinformationassignment.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(contigsassigned),sep="\t")



cat(paste0("Finished individual ",NAMES,"\n"))

