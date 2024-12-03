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


parser <- ArgumentParser(description= 'Gather blastn false positive results')
parser$add_argument('--inputRenv', '-I', help= 'Input R environment from summary99 summarise results rule')
parser$add_argument('--inputblastn_results', '-D', help= 'Input blastn_false_positive_check_results')
parser$add_argument('--outputpath', '-E', help= 'Output directory path, yes its the same as in input for this to')
parser$add_argument('--programdir', '-F', help= 'Program directory')
parser$add_argument('--inputpath', '-G', help= 'Input directory path')
parser$add_argument('--samplename', '-S', help= 'Sample name')



xargs2<- parser$parse_args()

# define base parameters and parameter variables for which analyses were run
NAMES <- xargs2$samplename
basepath <- xargs2$programdir
resultspath <- xargs2$outputpath
outtablespath <- paste0(basepath,resultspath)



blastnfalsepospresence <-readLines(xargs2$inputblastn_results)

if (length(blastnfalsepospresence) >=1) {
  Blastnfalseposhits <- read.table(file=xargs2$inputblastn_results, sep="\t",header=TRUE,row.names=NULL, fill=TRUE,quote="")
  paste0(NAMES," dim Blastnhits ", dim(Blastnfalseposhits))
}
if (length(blastnfalsepospresence) ==0) {
  Blastnfalseposhits <- as.data.frame(matrix(nrow=0,ncol=19))
  paste0(NAMES," Blastnhits returned no findings")
}

colnames(Blastnfalseposhits) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore","staxids", "stitle", "qcovhsp", "multiplesp","staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species","subspecies")

# This needs to be below the definition of the parser arguments pushed from this rule otherwise the parser args from gather_results_env just override these ones. Also saved
# args as xargs2 which should fix the issue to
Renv <- xargs2$inputRenv
load(Renv)


# I need this to change for the tables and figures for 99 compile
# allassignedfreqs


# Set as No, but set it as yes if the contigs are found 

allassignedfreqs <- allassignedfreqspreblastnfpcheckNas

allassignedfreqs$blastn_false_positive_check <- "No"
allassignedfreqs$blastn_alternate_superkingdom_id <- "NA"
allassignedfreqs$blastn_alternate_species <- "NA"
allassignedfreqs$blastn_alternate_subspecies <- "NA"
allassignedfreqs$blastn_alternate_percentident <- "NA"
allassignedfreqs$blastn_alternate_alignment_length <- "NA"


if (nrow(Blastnfalseposhits) >=1) {
  
  for (i in c(1:nrow(Blastnfalseposhits))) {
    
    
    grep(pattern = paste0(Blastnfalseposhits$qseqid[i],"$"),x=allassignedfreqs$contig) -> idx
    
    
    if (length(idx) >=1) {
      #print(paste0(" Match observed in ", Blastnfalseposhits$qseqid[i], allassignedfreqs$contig[idx], " index location ", idx, "species " , allassignedfreqs$subspecies[idx]))
      allassignedfreqs$blastn_false_positive_check[idx] <- "Yes"
      allassignedfreqs$blastn_alternate_superkingdom_id[idx] <- Blastnfalseposhits$superkingdom[i]
      allassignedfreqs$blastn_alternate_species[idx] <- Blastnfalseposhits$species[i]
      allassignedfreqs$blastn_alternate_subspecies[idx] <- Blastnfalseposhits$subspecies[i]
      allassignedfreqs$blastn_alternate_percentident[idx] <- Blastnfalseposhits$pident[i]
      allassignedfreqs$blastn_alternate_alignment_length[idx] <- Blastnfalseposhits$length[i]
      
    }
    
    
    
  }
  
}

allassignedfreqs2 <- subset(allassignedfreqs, !(is.na(allassignedfreqs$contigassignment) & (allassignedfreqs$blastn_alternate_superkingdom_id=="NA")))
allassignedfreqs <- allassignedfreqs2 





# I also want to output a new version for the output 20 and 100 viruses 
# Viraltop100


write.table(allassignedfreqs,file=(paste0(outtablespath,NAMES,"_summarycontighits_assigned_assembly_including_blastn_false_positive_check.txt")),sep="\t",row.names=FALSE,quote = FALSE)



# Need to prep the viral species counts to 
# Start with freq summary results

freqsummarynona$blastn_false_positive_check <- "No"
freqsummarynona$blastn_alternate_superkingdom_id <- "NA"
freqsummarynona$blastn_alternate_species <- "NA"
freqsummarynona$blastn_alternate_subspecies <- "NA"
freqsummarynona$blastn_alternate_percentident <- "NA"
freqsummarynona$blastn_alternate_alignment_length <- "NA"


if (nrow(Blastnfalseposhits ) >=1) {
  for (i in c(1:nrow(Blastnfalseposhits))) {
    
    
    grep(pattern = paste0(Blastnfalseposhits$qseqid[i],"$"),x=freqsummarynona$contig) -> idx
    
    
    if (length(idx) >=1) {
      
      freqsummarynona$blastn_false_positive_check[idx] <- "Yes"
      freqsummarynona$blastn_alternate_superkingdom_id[idx] <- Blastnfalseposhits$superkingdom[i]
      freqsummarynona$blastn_alternate_species[idx] <- Blastnfalseposhits$species[i]
      freqsummarynona$blastn_alternate_subspecies[idx] <- Blastnfalseposhits$subspecies[i]
      freqsummarynona$blastn_alternate_percentident[idx] <- Blastnfalseposhits$pident[i]
      freqsummarynona$blastn_alternate_alignment_length[idx] <- Blastnfalseposhits$length[i]
      
    }
    
    
    
  }
}

contigsspecies$contigs_assigned_to_species <- 0
contigsspecies$false_positive_blastn_test_undertaken <- "No"
contigsspecies$alternately_assigned_contigs <- "0"
contigsspecies$top_alternate_assigned_superkingdom <- "NA"
contigsspecies$top_alternate_assigned_species <- "NA"
contigsspecies$top_alternate_assigned_subspecies <- "NA"
contigsspecies$alternate_assigned_average_percent_ident <- "NA"
contigsspecies$alternate_assigned_average_alignment_length <- "NA"



if (nrow(Blastnfalseposhits ) >=1) {
  
  
  for (i in c(1:length(species_idvec))) {
    
    contigssubset<- subset(freqsummarynona,freqsummarynona$species == species_idvec[i])
    
    contigsspecies[i,2:9] <- contigssubset[1,4:11]
    contigsspecies[i,1] <- sum(contigssubset$freq)
    
    if (!is.na(contigssubset$percentident[1])) {
      contigsspecies[i,10] <- mean(contigssubset$percentident)
      contigsspecies[i,11] <- min(contigssubset$percentident)
      contigsspecies[i,12] <- max(contigssubset$percentident)
      contigsspecies[i,13] <- mean(contigssubset$contigalignlength)
      contigsspecies[i,14] <- nrow(contigssubset)
    }
    
    if (sum(contigssubset$blastn_false_positive_check == "Yes")>=1){
      
      contigsspecies[i,15] <- "Yes"
      contigsspecies[i,16] <-(sum(contigssubset$blastn_false_positive_check == "Yes"))
      
      
      counts <- contigssubset %>%
        filter(!(blastn_alternate_superkingdom_id=="NA")) %>%
        count(blastn_alternate_superkingdom_id, sort = TRUE)
      top_countsuperkingdom <- counts[1,1]
      
      
      
      counts <- contigssubset %>%
        filter(!(blastn_alternate_species=="NA")) %>%
        count(blastn_alternate_species, sort = TRUE)
      top_countspecies <- counts[1,1]
      
      counts <- contigssubset %>%
        filter(!(blastn_alternate_subspecies=="NA")) %>%
        count(blastn_alternate_subspecies, sort = TRUE)
      top_countsubspecies <- counts[1,1]
      
      
      
      avg_idents <- contigssubset %>%
        filter(!(blastn_alternate_percentident=="NA"))
      
      avg_identity_alt <- mean(as.numeric(avg_idents$blastn_alternate_percentident))
      
      
      
      
      avg_lengths <- contigssubset %>%
        filter(!(blastn_alternate_alignment_length=="NA"))
      
      avg_length_alt <- mean(as.numeric(avg_lengths$blastn_alternate_alignment_length))
      
      
      
      contigsspecies[i,17] <- top_countsuperkingdom
      contigsspecies[i,18] <- top_countspecies
      contigsspecies[i,19] <- top_countsubspecies
      contigsspecies[i,20] <- avg_identity_alt
      contigsspecies[i,21] <- avg_length_alt
      
      
      
    }
    
  }
  
  
}

Viralspecies <- subset(contigsspecies,contigsspecies$superkingdom=="Viruses")



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





write.table(Viraltop100,file=(paste0(outtablespath,NAMES,"_top100Viralhits_contigs_included_false_positive_check.txt")),sep="\t",row.names=FALSE,quote = FALSE)
write.table(Viraltop20,file=(paste0(outtablespath,NAMES,"_top20Viralhits_contigs_included_false_positive_check.txt")),sep="\t",row.names=FALSE,quote = FALSE)
