# Script for subsetting and extracting valuble information from Diamond blast
# Based off the script using in step 8 for extracting hits from contigs but will first calculate dram results and then add inferred contigs from blastx
# Diamond blast returns the top X (usually set to 10 to 20) hits for each contig. 
# I want to use this information to assign contigs with a putative assignment to species
# or if not possible, to the most recent ancestral grouping
# A quick method for this is to just go by what the top hit is and compile. This has some issues though
# If i took the top hit and the top hit has a vague tax id (e.g., to genus or family predicted) 
# Then that is little information, it might also be that the e values between them are practically
# identical (e.g.,  what is the functional difference between 1.5*10^-25 and 1.3*10^25)
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
library(phylotools)

parser <- ArgumentParser(description= 'Informing Diamond blast')

parser$add_argument('--inputdiamond', '-i', help= 'I am the input diamond file')
parser$add_argument('--output', '-o', help= 'Output file for summary of genomad returns')
parser$add_argument('--programdir ', '-p', help= 'working program directory')
parser$add_argument('--savdir ', '-s', help= 'working directory for saving contig matches')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--Acc', '-A', help= 'Accessiontaxa filepath')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')
parser$add_argument('--Log', '-l', help= 'Log of data')
parser$add_argument('--FDRrates', '-F', help= 'Estimated False discovery rates for genomad hits')

xargs<- parser$parse_args()

Diamond_output <- read.table(xargs$inputdiamond, header = FALSE, sep = "\t",fill=TRUE, quote="")

outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$savdir
outtablespath <- paste0(outtablespath1,outtablespath2)


FDRrates <- read.table(xargs$FDRrates, header = FALSE, sep = "\t",fill=TRUE)




NAMES <- xargs$name
Accessionfile <- xargs$Acc
AccessionNamenode <- xargs$Accnode
assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
assigned_contigs <- paste0(assigned_contigs1,assigned_contigs2)



# give proper column names
colnames(Diamond_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","qstart", "qend" ,"sstart", "send", "staxids", "stitle", "qcovhsp")

cat(paste0("Succesfully read in Diamond file for sample", NAMES,"\n"))

cat("Head 10 \n ")
print(Diamond_output[1:10,])


# This section checks to see whether the taxid of the optimal contig hits is one or multiple taxids
# This is important for downstream knowning closest species 
# but also important for actually generating larger taxonomic structure in the next step

# 
# Cut find optimal hit from here
# Lines 103-222 init. Lines 10-128 in overflow

#
# Temp redirect assuming filter done (may become permanent if filter removed)
contigsassigned <- Diamond_output

contigsassigned[,14] <- NA

# rd 1 generate tax id's from accessions

contigsassigned[,14] <- foreach(i = 1:nrow(contigsassigned), .combine = 'rbind') %dopar% {
  
  contigsassigned[,14] <- taxonomizr::accessionToTaxa(contigsassigned$sseqid[i], sqlFile=Accessionfile, version = "version")
}

# rd 2. Many accessions don't work (old, outdated? maybe something else)
# generate taxids from name in the stitle

contigsassigned[,15] <- "name"

for (i in c(1:nrow(contigsassigned))) {
  

  string <- str_extract(string = contigsassigned$stitle[i],pattern = "\\[.*\\]")
  string2 <- gsub(x=string, pattern = "[",replacement = "",fixed = TRUE)
  string3 <- gsub(x=string2, pattern = "]",replacement = "",fixed = TRUE)
  string4 <- gsub(x=string3,pattern = " sp." ,replacement= "",fixed = TRUE)
  contigsassigned[i,15] <- string4
  
}

for ( i in c (1:nrow(contigsassigned))) {
  
  
  if (is.na(contigsassigned$V14[i])) {
    
    
    getId(contigsassigned$V15[i], sqlFile = AccessionNamenode, onlyScientific = FALSE)
    
  }
  
}


contigsassigned[,14] <- foreach(i = 1:nrow(contigsassigned), .combine = 'rbind') %dopar% {
  
  contigsassigned[,14] <- taxonomizr::accessionToTaxa(contigsassigned$sseqid[i], sqlFile=Accessionfile, version = "version")
}



contigsassigned[,16:22] <- NA


colnames(contigsassigned) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","qstart", "qend" ,"sstart", "send", "staxids", "stitle", "qcovhsp","taxid","name","Kingdom","phylum","class","order","family","genus","species")

contigsassigned[,16:22]  <- foreach(i = 1:nrow(contigsassigned), .combine = 'rbind') %dopar% {
  
  contigsassigned[i,16:22]  <- taxonomizr::getTaxonomy(contigsassigned$taxid[i], sqlFile=AccessionNamenode)
  
}

contigsassigned$FDRrate <- 0
contigsassigned$Viral_hallmarks <- 0

for (i in c(1:nrow(FDRrates))) {
  
  
  if (length(grep(FDRrates[i,1],contigsassigned$qseqid)) >=1 ) {
    
    idx <- grep(FDRrates[i,1],contigsassigned$qseqid)
    
    contigsassigned$FDRrate[idx] <- as.numeric(FDRrates[i,8])
    contigsassigned$Viral_hallmarks[idx] <- as.numeric(FDRrates[i,9])
  }
  
  
}




uniquetax <- unique(contigsassigned$taxid)


if (!is.null(nrow(uniquetax))) {
taxlist <- list()

for ( i in c(1:nrow(uniquetax))) {
  
  grep(uniquetax[i,],contigsassigned$taxid) -> idx
  
  sampledata <- contigsassigned[idx,]
  
  taxlist[[i]] <- sampledata
  
}


sum(contigsassigned$Viral_hallmarks >=1)



d <- Sys.time()

cat(paste0(NAMES," Finished taxid conversion to taxonomy ", "\n"))
cat(paste0(Sys.time(), "\n"))


length(taxlist)

# store lowest confidence (only 1 contig assigned, no contigs)
taxlist1 <- list()

# taxlist2 is going to store all low-mid confidence (more than 1 contig but found at the same/overlapping parts of only one gene but no hallmarks)
taxlist2 <- list()


# taxlist3 is going to store all mid confidence (more than 1 contig which are spread across different genes or distinct regions of one gene but no hallmarks)
taxlist3 <- list()

# taxlist4 is going to store highest confidence viruses (At least one contig and the presence of at least one viral hallmark or at least 2 genes of high pw identity to an existing virus  or at least one contig at very high pw identity)
taxlist4 <- list()

summaryvirushits <- as.data.frame(matrix(nrow=length(taxlist),ncol = 8))
colnames(summaryvirushits) <- c("Viral_species", "N_contigs_assigned", "prot_ids_aligned","unique_genes_aligned_to","average_identity","number_viral_hallmarks", "avg_FDR_rate","certainty_category")
h=1
j=1
k=1
l=1

for (i in c(1:length(taxlist))) {
  
  multimarker=FALSE
  assignement = "No"
  as.data.frame(taxlist[i]) -> virusX 
  if (!is.na(virusX$pident[1])) {
  ids <- gsub(pattern="_[0-9]$",replacement="",x= virusX$qseqid)
  ids <- unique(ids)
  # remove the accession code 
  genes <- gsub(pattern = "^.*\\.1 ",replacement ="" ,x = virusX$stitle)
  genes <- gsub(pattern = "^.*\\.2 ",replacement ="" ,x = genes)
  # remove any MAG code
  genes <- gsub(pattern = "MAG: ",replacement ="" ,x = genes)
  # Remove incorrectly added .
  genes <- gsub(pattern = "\\.",replacement ="" ,x = genes)
  # Remove incorrectly added double spacing
  genes <- gsub(pattern = "  ",replacement =" " ,x = genes)
  # Remove partialin case one contig hits to a partial protein and another hits to the same protein but from a source that was complete
  genes <- gsub(pattern = ", partial",replacement ="" ,x = genes)
  # as above but reverse
  genes <- gsub(pattern = ", complete",replacement ="" ,x = genes)
  # redo but assuming no ,
  genes <- gsub(pattern = " partial",replacement ="" ,x = genes)
  # redo but assuming no ,
  genes <- gsub(pattern = " complete",replacement ="" ,x = genes)
  # remove species name if also added to the gene tag
  genes <- gsub(pattern = " \\[.*$",replacement ="" ,x = genes)
  #number of unique genes present
  length(unique(genes))
  # extract a gene (doesn't matter which because hits with multiple genes aren't filtered by these values)
  grep(pattern = genes[1],x = virusX[,12], fixed=TRUE) -> geneidx 
  # subset to only rows by that gene
  virusXgene1 <- virusX[geneidx,]
  # Extract the min and max values of alignement to that gene 
  grep(paste0("^",min(virusXgene1$sstart),"$"),virusXgene1$sstart) -> minval
  grep(paste0("^",max(virusXgene1$send),"$"),virusXgene1$send) -> maxval
  

  # Check if they are in the same rows.
  virusXgene1minmax <- virusXgene1[c(minval,maxval),]

  ids2 <- gsub(pattern="_[0-9]$",replacement="",x= virusXgene1minmax$qseqid)
  ids2 <- unique(ids2)
  
  if ( (length(unique(genes)) >= 2)) {
  
    multimarker=TRUE
  }
    if ( (minval[1] != maxval[1] & length(ids2) >= 2) ) {
    # second largest end value
      
    #largest end value
    largeend <- virusXgene1$send[maxval[1]]
    largestart <- virusXgene1$sstart[maxval[1]]
    # smallest end values
    smallend <- virusXgene1$send[minval[1]]
    smallstart <- virusXgene1$sstart[minval[1]]
    
    if (largeend >= (smallend + 150) | smallstart <= (largestart - 150) ) {

      multimarker=TRUE
    }
  } 
  
  # what is the average pident. 
  avgpident <- mean(virusX$pident)
  contig1length <- virusX$length[1]
  
  # Are hallmarks present or is percent identity very high? If so, highest confidence assignment
  if ((sum(virusX$Viral_hallmarks) >= 2) | (sum(virusX$Viral_hallmarks) >= 1 & avgpident >=75) | ((length(unique(genes)) >=3 & avgpident >=85)) | (length(unique(ids)) >= 3 & avgpident >=90) | (avgpident >=95 & contig1length >=1000) ) {
    taxlist4[[j]] <- virusX
    j=j+1
    
    summaryvirushits$Viral_species[i] <- virusX$species[1]
    
    ids <- gsub(pattern="_[0-9]$",replacement="",x= virusX$qseqid)
    ids <- unique(ids)
    summaryvirushits$N_contigs_assigned[i] <- length(ids)
    summaryvirushits$prot_ids_aligned[i] <- length(unique(virusX$sseqid))
    
    summaryvirushits$unique_genes_aligned_to[i] <- length(unique(genes))
    summaryvirushits$average_identity[i] <- mean(virusX$pident)
    summaryvirushits$number_viral_hallmarks[i] <- sum(virusX$Viral_hallmarks)
    summaryvirushits$avg_FDR_rate[i] <- mean(virusX$FDRrate)
    summaryvirushits$certainty_category[i] <- "Most certain (4)"
    
    assignement = "assigned"
  }
  
  
  
  if ((sum(virusX$Viral_hallmarks) == 0) & (length(ids) >=3 | avgpident >=85) & (length(unique(genes)) >=3 | multimarker ==TRUE) & assignement != "assigned") {
  
    taxlist3[[k]] <- virusX
    k=k+1
    
    summaryvirushits$Viral_species[i] <- virusX$species[1]
    
    ids <- gsub(pattern="_[0-9]$",replacement="",x= virusX$qseqid)
    ids <- unique(ids)
    summaryvirushits$N_contigs_assigned[i] <- length(ids)
    summaryvirushits$prot_ids_aligned[i] <- length(unique(virusX$sseqid))
    
    
    summaryvirushits$unique_genes_aligned_to[i] <- length(unique(genes))
    summaryvirushits$average_identity[i] <- mean(virusX$pident)
    summaryvirushits$number_viral_hallmarks[i] <- sum(virusX$Viral_hallmarks)
    summaryvirushits$avg_FDR_rate[i] <- mean(virusX$FDRrate)
    summaryvirushits$certainty_category[i] <- "High confidence (3)"
    assignement = "assigned"
  
  }

  
  if ((sum(virusX$Viral_hallmarks) == 0) & (length(ids) >=2) & (multimarker ==FALSE) & assignement != "assigned") {
    
    taxlist2[[l]] <- virusX
    l=l+1
    
    summaryvirushits$Viral_species[i] <- virusX$species[1]
    
    ids <- gsub(pattern="_[0-9]$",replacement="",x= virusX$qseqid)
    ids <- unique(ids)
    summaryvirushits$N_contigs_assigned[i] <- length(ids)
    summaryvirushits$prot_ids_aligned[i] <- length(unique(virusX$sseqid))
    
    
    summaryvirushits$unique_genes_aligned_to[i] <- length(unique(genes))
    summaryvirushits$average_identity[i] <- mean(virusX$pident)
    summaryvirushits$number_viral_hallmarks[i] <- sum(virusX$Viral_hallmarks)
    summaryvirushits$avg_FDR_rate[i] <- mean(virusX$FDRrate)
    summaryvirushits$certainty_category[i] <- "Mid confidence (2)"
    
    assignement = "assigned"
  }
  
  

  if ((sum(virusX$Viral_hallmarks) == 0) & (length(ids) <=1) & assignement != "assigned") {
    
    taxlist1[[h]] <- virusX
    h=h+1
    
    summaryvirushits$Viral_species[i] <- virusX$species[1]
    
    ids <- gsub(pattern="_[0-9]$",replacement="",x= virusX$qseqid)
    ids <- unique(ids)
    summaryvirushits$N_contigs_assigned[i] <- length(ids)
    summaryvirushits$prot_ids_aligned[i] <- length(unique(virusX$sseqid))
 
    
    summaryvirushits$unique_genes_aligned_to[i] <- length(unique(genes))
    summaryvirushits$average_identity[i] <- mean(virusX$pident)
    summaryvirushits$number_viral_hallmarks[i] <- sum(virusX$Viral_hallmarks)
    summaryvirushits$avg_FDR_rate[i] <- mean(virusX$FDRrate)
    summaryvirushits$certainty_category[i] <- "Low confidence (1)"
    
    assignement = "assigned"
  }
  
  
  
  
}
  
}


# low quality
taxlist1

#low-mid quality
taxlist2


# mid quality
taxlist3

# High quality
taxlist4


# Will need to write the summary table. 
# will need to write each individual viral species for high and highest confidence
# Will cat low and mid and export as other. 


if (length(taxlist3) >=1) {
for ( i in (c(1:length(taxlist3)))) {
  
  taxlist3[[i]] -> viralsptable
  
  viralsptable$species[1] -> spname
  spname <- gsub(" ","_",spname)
  spname <- gsub("/","_",spname)
  
  exportname <- paste0("/",NAMES,"_Viral_contigs_for_", spname)
  
  write.table(viralsptable,file = paste0(outtablespath,exportname,".tsv"),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)

  }

}

if (length(taxlist4) >=1) {

for ( i in (c(1:length(taxlist4)))) {
  
  taxlist4[[i]] -> viralsptable
  
  viralsptable$species[1] -> spname
  spname <- gsub(" ","_",spname)
  spname <- gsub("/","_",spname)
  
  exportname <- paste0("/",NAMES,"_Viral_contigs_for_", spname)
  
  write.table(viralsptable,file = paste0(outtablespath,exportname,".tsv"),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)
  
}
}



tax1values <- bind_rows(taxlist1)
tax2values <- bind_rows(taxlist2)
tax3values <- bind_rows(taxlist3)
tax4values <- bind_rows(taxlist4)



exportname <- paste0("/",NAMES,"_uncertain_virus_calls_all")
exporthighest <- paste0("/",NAMES,"_most_certain_virus_calls_all")
exporthigh <- paste0("/",NAMES,"_high_and_certain_virus_calls_all")
uncertaintaxa <- rbind(tax1values,tax2values)
highandcertaincertaintaxa <- rbind(tax3values,tax4values)

# low confidence
write.table(uncertaintaxa,file = paste0(outtablespath,exportname,".tsv"),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)

# High confidence
write.table(tax4values,file = paste0(outtablespath,exporthighest,".tsv"),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(highandcertaincertaintaxa,file = paste0(outtablespath,exporthigh,".tsv"),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)



summaryvirushits <- summaryvirushits[complete.cases(summaryvirushits),]

write.table(summaryvirushits,file = paste0(assigned_contigs),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)

}


if (is.null(nrow(uniquetax))) {
summaryvirushits <- as.data.frame(matrix(nrow=2, ncol=2))

write.table(summaryvirushits,file = paste0(assigned_contigs),quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)


}