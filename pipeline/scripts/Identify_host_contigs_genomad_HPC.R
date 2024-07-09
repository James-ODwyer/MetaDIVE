# Investigate potential of including host contigs assigned from summary in downstream genomad analysis (save on another diamond)

# Looking at outputs, it won't work. I need the multiple hits per tag to see if it's widespread. 

# Try the contigs output from blastn 


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
parser$add_argument('--output', '-o', help= 'Output file for host contigs')
parser$add_argument('--programdir ', '-p', help= 'working program directory')
parser$add_argument('--savdir ', '-s', help= 'working directory for saving contig matches')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--Acc', '-A', help= 'Accessiontaxa filepath')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')
parser$add_argument('--Log', '-l', help= 'Log of data')
parser$add_argument('--hostsp', '-C', help= 'host species if determined')
parser$add_argument('--contighostsp', '-P', help= 'host species determined through contig blastx')

xargs<- parser$parse_args()

Diamond_output <- read.table(xargs$inputdiamond, header = FALSE, sep = "\t",fill=TRUE)

outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$savdir
outtablespath <- paste0(outtablespath1,outtablespath2)




hostspeciesLines <- readLines(xargs$hostsp)

if (length(hostspeciesLines) >= 1 ) {
  hostspeciescheck <- "Present"
  hostspecies <- read.table(xargs$hostsp)
  
  
}


if (length(hostspeciesLines) < 1 ) {
  
  hostspeciescheck <- "Absent"
  
  
}

hostspeciescontigspth <- xargs$contighostsp





NAMES <- xargs$name
Accessionfile <- xargs$Acc
AccessionNamenode <- xargs$Accnode
assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
assigned_contigs <- paste0(assigned_contigs1,assigned_contigs2)


log <- xargs$Log


hostspeciescontigs<- read.table(hostspeciescontigspth)


#This part I need to check manually before I can succesfully incorporate it into the pipeline. Will it import as quote, will it import from table? how will it import? etc etc
# If it all works, I can skip the assigning phase and just jump straight to removing the host contigs based in the taxid instead. 


# Input finished Diamond (when slotting into Diamond R script, can ignore)

colnames(Diamond_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")


contigsassigned <- Diamond_output






# So we have a completed Diamond assigned andtax id listed blast.

# Now call in the host species to identify family level.

# If host was chosen to be detected, this is easy

  
  
  # Identify what are probable host using identified host data
  
  
  # Use host id genus to collect family from name without needing to generate tax info again (can also use species it doesn't matter)

  # May need to change back to original species id of I have removed everything with revised host contigs step and don't generate any species data from contigs now



grep(hostspeciescontigs$V1,contigsassigned$staxids) -> idx
  

contigsassignedhostonly <-   contigsassigned[idx,]
contigsassignednohost <-   (contigsassigned[-idx,])
  

contigsassignedhostonly$qseqid -> hostcontigs
# If not we will need to infer likely host from contigs returned here. 
# This will be done as the top returned family. It will as far as I can see always be the host family given expected DNA concentrations for the host
# Although now I have prefiltered out host contigs in earlier steps maybe the ratios will be more even. Will have to inspect manually.


# Now per Diamond hit if over 80% homology is found to 2+ distict proteins of the family. the contig is considered mosquito and will be removed
# before genomad


write.table(hostcontigs,file=assigned_contigs,sep="\t", row.names=FALSE, quote=FALSE,col.names=FALSE)



