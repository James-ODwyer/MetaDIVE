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

# --- standardize to 25-col schema (adds alternate_* columns if missing) ---
contig_cols25 <- c(
  "qseqid","sseqid","pident","length","evalue","bitscore",
  "staxids","stitle","qcovhsp","multiplesp","staxidreduced",
  "superkingdom","phylum","class","order","family","genus","species","subspecies",
  "alternate_genus1","alternate_species1",
  "alternate_genus2","alternate_species2",
  "alternate_genus3","alternate_species3"
)

std_to_25 <- function(df) {
  miss <- setdiff(contig_cols25, names(df))
  for (m in miss) df[[m]] <- NA
  df <- df[, contig_cols25, drop = FALSE]
  df
}

Blastnfalseposhits <- std_to_25(Blastnfalseposhits)

# Coerce numeric fields we use later
suppressWarnings({
  Blastnfalseposhits$pident  <- as.numeric(Blastnfalseposhits$pident)
  Blastnfalseposhits$length  <- as.numeric(Blastnfalseposhits$length)
  Blastnfalseposhits$qcovhsp <- as.numeric(Blastnfalseposhits$qcovhsp)
})
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

# Reorder the df so that the data has the secondary hits last!!! Still TODO
allassignedfreqs


# I also want to output a new version for the output 20 and 100 viruses 
# Viraltop100


write.table(allassignedfreqs,file=(paste0(outtablespath,NAMES,"_summarycontighits_assigned_assembly_including_blastn_false_positive_check.txt")),sep="\t",row.names=FALSE,quote = FALSE)



# Need to prep the viral species counts to 
# Start with freq summary results

alt_cols <- c("alternate_genus1","alternate_species1",
              "alternate_genus2","alternate_species2",
              "alternate_genus3","alternate_species3")
for (cc in alt_cols) {
  if (!cc %in% names(freqsummarynona)) freqsummarynona[[cc]] <- NA
}


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

alt_cols_cs <- c("alternate_genus1","alternate_species1",
                 "alternate_genus2","alternate_species2",
                 "alternate_genus3","alternate_species3")
for (cc in alt_cols_cs) {
  if (!cc %in% names(contigsspecies)) contigsspecies[[cc]] <- NA
}

# most frequent non-NA/non-blank value, default to "NONE"
mode_nonempty <- function(v) {
  v <- v[!is.na(v) & v != ""]
  if (length(v) == 0) return("NONE")
  names(sort(table(v), decreasing = TRUE))[1]
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
      # write to named base columns
      contigsspecies[i, "average_percent_ident"] <- mean(contigssubset$percentident)
      contigsspecies[i, "min_percent_ident"]     <- min(contigssubset$percentident)
      contigsspecies[i, "max_percent_ident"]     <- max(contigssubset$percentident)
      contigsspecies[i, "length"]                <- mean(contigssubset$contigalignlength)
    }
    
    # how many contigs contributed to this species
    contigsspecies[i, "contigs_assigned_to_species"] <- nrow(contigssubset)
    
    # fill alternates for this species by most frequent across its contigs
    contigsspecies[i, "alternate_genus1"]   <- mode_nonempty(contigssubset$alternate_genus1)
    contigsspecies[i, "alternate_species1"] <- mode_nonempty(contigssubset$alternate_species1)
    contigsspecies[i, "alternate_genus2"]   <- mode_nonempty(contigssubset$alternate_genus2)
    contigsspecies[i, "alternate_species2"] <- mode_nonempty(contigssubset$alternate_species2)
    contigsspecies[i, "alternate_genus3"]   <- mode_nonempty(contigssubset$alternate_genus3)
    contigsspecies[i, "alternate_species3"] <- mode_nonempty(contigssubset$alternate_species3)
    
    # summarize Blastn FP alternates
    if (sum(contigssubset$blastn_false_positive_check == "Yes") >= 1) {
      contigsspecies[i, "false_positive_blastn_test_undertaken"] <- "Yes"
      contigsspecies[i, "alternately_assigned_contigs"] <- sum(contigssubset$blastn_false_positive_check == "Yes")
      
      # top alternate superkingdom
      counts <- contigssubset %>%
        filter(!(blastn_alternate_superkingdom_id == "NA")) %>%
        count(blastn_alternate_superkingdom_id, sort = TRUE)
      top_countsuperkingdom <- if (nrow(counts) >= 1) counts$blastn_alternate_superkingdom_id[1] else "NA"
      
      # top alternate species
      counts <- contigssubset %>%
        filter(!(blastn_alternate_species == "NA")) %>%
        count(blastn_alternate_species, sort = TRUE)
      top_countspecies <- if (nrow(counts) >= 1) counts$blastn_alternate_species[1] else "NA"
      
      # top alternate subspecies
      counts <- contigssubset %>%
        filter(!(blastn_alternate_subspecies == "NA")) %>%
        count(blastn_alternate_subspecies, sort = TRUE)
      top_countsubspecies <- if (nrow(counts) >= 1) counts$blastn_alternate_subspecies[1] else "NA"
      
      # averages for FP alternates
      avg_identity_alt <- contigssubset %>%
        filter(!(blastn_alternate_percentident == "NA")) %>%
        mutate(bapi = as.numeric(blastn_alternate_percentident)) %>%
        { if (nrow(.) > 0) mean(.$bapi, na.rm = TRUE) else NA_real_ }
      
      avg_length_alt <- contigssubset %>%
        filter(!(blastn_alternate_alignment_length == "NA")) %>%
        mutate(bal = as.numeric(blastn_alternate_alignment_length)) %>%
        { if (nrow(.) > 0) mean(.$bal, na.rm = TRUE) else NA_real_ }
      
      contigsspecies[i, "top_alternate_assigned_superkingdom"]        <- top_countsuperkingdom
      contigsspecies[i, "top_alternate_assigned_species"]             <- top_countspecies
      contigsspecies[i, "top_alternate_assigned_subspecies"]          <- top_countsubspecies
      contigsspecies[i, "alternate_assigned_average_percent_ident"]   <- avg_identity_alt
      contigsspecies[i, "alternate_assigned_average_alignment_length"]<- avg_length_alt
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
