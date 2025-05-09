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
library("readr")
library(future.apply)

#args <- commandArgs(TRUE)
parser <- ArgumentParser(description= 'Informing Diamond blast')

parser$add_argument('--inputdiamond', '-i', help= 'I am the input diamond file')
parser$add_argument('--inputvirusall', '-b', help= 'I am the csv results for virus_all')
parser$add_argument('--output', '-o', help= 'Output file for combined csvs')
parser$add_argument('--programdir ', '-p', help= 'working program directory')
parser$add_argument('--savdir ', '-s', help= 'working directory for saving contig matches')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--threads', '-t', help= 'Number of threads')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')


xargs<- parser$parse_args()



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


# Apply the function to all lines

resultspathcompile <- xargs$savdir

modified_lines <- lapply(Diamond_lines, replace_newline_within_quotes)
modified_lines  <- unlist(modified_lines)
modified_lines <- modified_lines[sapply(modified_lines, function(x) length(strsplit(x, "\t")[[1]]) > 1)]
# Read the modified lines using read.table
Diamond_output <- read.table(text = modified_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")

# Replace the unique character sequence back to '\n' within the dataframe
Diamond_output[] <- lapply(Diamond_output, function(col) gsub("###NEWLINE###", "\n", col))

rm(modified_lines,Diamond_lines)


outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$abundances
outtablespath <- paste0(outtablespath1,outtablespath2)
AccessionNamenode <- xargs$Accnode

NAMES <- xargs$name
ncores <- xargs$threads


assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
assigned_contigs <- paste0(assigned_contigs1,assigned_contigs2)

contigsvdir <-xargs$savdir
svcontigname <-xargs$savcontig

log <- xargs$Log

plan(multisession, workers = ncores)

# Set to correct directory.
# sink all cats and prints to the snakemake log files 
#sink(log)


if (!file.exists(outtablespath)){
  
  dir.create(outtablespath, recursive = FALSE, mode = "0777")
}


# give proper column names
colnames(Diamond_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")

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
# Later update here 30/05/2024, Diamond appears to fail on proteins specifically PDB. Have implemented
# a rapid taxid check using TaxonKit as a bash script in the rule. Problem should be gone
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


if ( nrow(uniquespmissing) >=1 ) {
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
process_row <- function(i) {
  row <- list()
  # assign taxid to row variable
  row[1] <- contigsassignedunique$staxidreduced[i]
  # Get standard taxonomy
  row[2:8] <- taxonomizr::getTaxonomy(contigsassignedunique$staxidreduced[i], sqlFile = AccessionNamenode)
  
  # Check for subspecies allocation
  values <- taxonomizr::getRawTaxonomy(contigsassignedunique$staxidreduced[i], sqlFile = AccessionNamenode)
  row[9] <- if (!is.null(values[[1]][1])) values[[1]][1] else row[8]
  
  # Handle missing taxid
  if (is.na(row[2])) {
    spname <- str_extract_all(contigsassignedunique$stitle[i], "\\[(.*?)\\]")[[1]]
    spname2 <- str_replace_all(spname, "\\[|\\]", "")
    value <- taxonomizr::getId(spname2, sqlFile = AccessionNamenode)
    
    row[2:8] <- taxonomizr::getTaxonomy(value, sqlFile = AccessionNamenode)
    
    values <- taxonomizr::getRawTaxonomy(value, sqlFile = AccessionNamenode)
    row[9] <- if (!is.null(values[[1]][1])) values[[1]][1] else row[8]
  }
  
  return(unlist(row))
}

# Apply the Function in Parallel
taxidsunique <- future_lapply(1:nrow(contigsassignedunique), process_row)

# Combine Results into Data Frame
taxidsunique <- do.call(rbind, taxidsunique)

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

rm(Diamond_output)

cat(paste0(NAMES," Finished reindexing of taxa results ", "\n"))
cat(paste0(Sys.time(), "\n"))

contigsassigned$stitle = substr(contigsassigned$stitle,1,50)
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
Eukaryotes <- subset(contigsassigned, contigsassigned$superkingdom=="Eukaryota")
Bacteria <- subset(contigsassigned, contigsassigned$superkingdom=="Bacteria")
Viruses <- subset(contigsassigned, contigsassigned$superkingdom=="Viruses")


totalcontigsassigned <- nrow(contigsassigned)

Eukaryotessinglesp <- Eukaryotes
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


Bacteriasinglesp <- Bacteria
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


virus_all <- read.csv(file = xargs$inputvirusall, header = TRUE, stringsAsFactors = FALSE)



virus_all$additional_diverged_reads <- 0
virus_all$total_reads_assigned <- as.numeric(virus_all$total_reads_assigned)

Virusessinglespspcounts$Frequency <- as.numeric(Virusessinglespspcounts$Frequency)

Virusessinglespspcounts$Species <- as.character(Virusessinglespspcounts$Species)

# Perform a left join to bring in Frequency where species matches
virus_all <- virus_all %>%
  left_join(
    Virusessinglespspcounts %>% select(Species, Frequency), 
    by = c("species" = "Species")
  ) %>%
  mutate(additional_diverged_reads = coalesce(Frequency, additional_diverged_reads)) %>%  # Replace only if a match is found
  select(-Frequency)  # Remove the extra Frequency column from the join


virus_all <- virus_all %>%
  mutate(
    total_reads_combined_diverged = coalesce(additional_diverged_reads, 0) + coalesce(total_reads_assigned, 0)
  )


  combined_virus_reads <- as.data.frame(matrix(nrow=0,ncol=2))
  colnames(combined_virus_reads) <- c("Reads","Species")
  combined_virus_contigs <- as.data.frame(matrix(nrow=0,ncol=2))
  colnames(combined_virus_contigs) <- c("Contigs","Species")
  

for (i in c(1:nrow(virus_all))) {
    
   
    workingcontigs <- subset(Viruses,Viruses$species==virus_all$species[i])
    
    contigsdf <- as.data.frame(matrix(nrow=nrow(workingcontigs),ncol=2))
    contigsdf$V1 <- workingcontigs[,1]
    
    if (nrow(contigsdf) >=1){
      contigsdf$V2 <- virus_all$species[i]
    }
    colnames(contigsdf) <- c("Contigs","Species")
    
    patterns <- unique(workingcontigs$contig)
    
    if (length(patterns)>=1) {
      pattern_regex <- paste(patterns, collapse = "|")
      grep(pattern=pattern_regex,x=reads_contigs_virus$contig) ->idx2
      
    }
    if (length(patterns)<1) {
      idx2 <- NULL 
    }
     
 
  combined_virus_contigs <- rbind(combined_virus_contigs,contigsdf)
    
  }

  combined_virus_contigs_unique <- unique(combined_virus_contigs)

  combined_virus_contigs_unique$Species <- gsub(" ", "_", combined_virus_contigs_unique$Species)

  directory <- resultspathcompile

  split_dfcontigs <- split(combined_virus_contigs_unique, combined_virus_contigs_unique$Species)
  names(split_dfcontigs) <- gsub("/", "_", names(split_dfcontigs))
  
  # Save each dataframe in the list as a separate text file
  for (category_value in names(split_dfcontigs)) {
    file_name <- paste0(directory, category_value, "/",category_value,"_divergent_reads_and_contigs.txt")
    contignames <- split_dfcontigs[[category_value]][, 1]
    write.table(split_dfcontigs[[category_value]], file = file_name, sep = "\t", row.names = FALSE,col.names=FALSE,quote = FALSE)
  }


  write.csv(virus_all, file=paste0(resultspathcompile,NAMES,"_virusall_sums_incl_diverged_reads.csv"), row.names = FALSE)

