library(dada2)
library(Biostrings)
library(stringr)
library(plyr)
library(dplyr)
library(argparse)
library(DT)
library(htmlTable)
library(htmltools)
library(rmarkdown)

parser <- ArgumentParser(description= 'Collecting_all_viral_results')

parser$add_argument('--inputvirusdetails', '-i', help= 'I am viral complexity summary text file ')
parser$add_argument('--name', '-N', help= 'Name of sample')
parser$add_argument('--threads', '-p', help= 'Number of threads (currently only one, not parallelised')
parser$add_argument('--contigs', '-C', help= 'contigs_file_generated')
parser$add_argument('--programdir', '-t', help= 'directory of analysis run (same as program_dir in config file)')
parser$add_argument('--savdir', '-A', help= 'directory path of where to save viral raw reads output table')
parser$add_argument('--Rdatas', '-r', help= ' I am the results Rdata environment from 99summary ')
parser$add_argument('--dofalseposcheck', '-f', help= ' Whether to subset to only contigs which passed blastn analysis ')
parser$add_argument('--allrawsassign', '-g', help= 'the complete blastn assignment table for raw reads')


xargs<- parser$parse_args()

viralcomplexity <- xargs$inputvirusdetails
NAMES <- xargs$name
basepath <- xargs$programdir
resultspath <- xargs$savdir
resultspathcompile <- xargs$savdir
contigspathing <- xargs$contigs
rdataenv <- xargs$Rdatas
threads <- xargs$threads
rawsassign <- xargs$allrawsassign
dofalsepossubset <- xargs$dofalseposcheck
outtablespathcompile <- paste0(basepath,resultspath)
outtablespath <- paste0(basepath,resultspath)



# working script for summary final virus reads detection

# The intention of this R script is to 
# 1. Generate complexity score for the contigs 
# 2. collate the raw reads and the contigs viral hits
# 3. Generate predicted certainty around the scores for certainty of occuring. 
# 4. Prepare for extracting all read/contig names for each virus to then be turned into Fastq files for storage in new folder 099_viral_reads

# Files will have column names still so file size is 434 (~12 columns with titles)
# I have set the complexity min to 500 wwhich should cover a bit of safety in case the names are changed a bit

load(rdataenv)



if (file.info(viralcomplexity)$size >= 330) {
  
	viral_raws_complexity_summary <- read.table(viralcomplexity, header=TRUE)
	rawsassigntable <- read.table(rawsassign, header=FALSE, sep="\t")

  
}

if (file.info(viralcomplexity)$size < 330) {
  
  viral_raws_complexity_summary <- as.data.frame(matrix(nrow=0, ncol=ncol(Combined_assigned_contigs)))
  colnames(viral_raws_complexity_summary) <- colnames(Combined_assigned_contigs)
	rawsassigntable <- read.table(rawsassign, header=FALSE,sep="\t")
}




# 1 Create complexity results for contigs
# For this I will need the fasta sequences from trinity/Megahit. 

# analysis may fail here! 
# Trinity returns fastas with a split line at 70 while megahit returns a wrapped line within the \n. How this will affect sequences remains to be seen 


sequences <- getSequences(contigspathing)
sequencenames <- names(sequences)
sequencenames <- sub(" .*", "", sequencenames)


complexity_scores <- seqComplexity(sequences, kmerSize=3)





seqs_plus_scores <- as.data.frame(cbind(sequencenames,complexity_scores))

Combined_assigned_contigs_viruses_only <- subset(Combined_assigned_contigs,Combined_assigned_contigs$superkingdom=="Viruses")

# If to replace viruses only if only false pos subsets are required

if ( dofalsepossubset == "confirmed") {
  
  Combined_assigned_contigs_viruses_only <- subset(Combined_assigned_contigs,Combined_assigned_contigs$superkingdom=="Viruses")
  
}

Combined_assigned_contigs_viruses_only$qseqid <- sub(" .*", "", Combined_assigned_contigs_viruses_only$qseqid)

seqscores_viruses <- seqs_plus_scores[seqs_plus_scores$sequencenames %in% Combined_assigned_contigs_viruses_only$qseqid, ]


# Start of large if loop for checking of contigs are empty

if(nrow(Combined_assigned_contigs_viruses_only) >=1) {
  Combined_assigned_contigs_viruses_only$complexity_score <- NA
  
  
  for (i in c(1:nrow(seqscores_viruses))) {
    
	grep(pattern=paste0("^", seqscores_viruses$sequencenames[i], "$"),x=Combined_assigned_contigs_viruses_only$qseqid) -> idx
    
    
	Combined_assigned_contigs_viruses_only$complexity_score[idx] <- seqscores_viruses$complexity_scores[i]
    
  }
  

  
  Combined_assigned_contigs_viruses_species_summary <- Combined_assigned_contigs_viruses_only %>%
    group_by(species) %>%
    summarize(count = n())
  
  
  Viraltop100 <- Viraltop100 %>%
    rename(reads_assigned_through_contigs = Frequency)
  
  Viraltop100$number_contigs_assigned <- NA
  Viraltop100$mean_complexity_of_contigs <- NA
  Viraltop100$contigs_assigned <- "yes"
  
  Combined_assigned_contigs_viruses_only$complexity_score <- as.numeric(Combined_assigned_contigs_viruses_only$complexity_score)



  
  for (i in c(1:nrow(Viraltop100))) {
    
    
    grep(pattern=Viraltop100$species[i],x=Combined_assigned_contigs_viruses_only$species) -> idx
    grep(pattern=Viraltop100$species[i],x=Combined_assigned_contigs_viruses_species_summary$species) -> idx2
    
    subsetspecies <- Combined_assigned_contigs_viruses_only[idx,]
    
    Viraltop100$number_contigs_assigned[i] <- Combined_assigned_contigs_viruses_species_summary$count[idx2]
    Viraltop100$mean_complexity_of_contigs[i] <- mean(subsetspecies$complexity_score)
    
    
    
  }
  
  Viraltop100$additional_raw_reads_assigned <- "no"
  Viraltop100$number_of_raw_reads_assigned <- 0
  Viraltop100$mean_complexity_of_raw_reads <- 0
  Viraltop100$mean_identity_of_raw_reads <- 0
  Viraltop100$mean_aligned_length_of_raw_reads <- 0
  

  # generate the complete viral data for the raws info added to the contigs. Do another grep loop
  # then subset the raws file to remove any species that are already in the contigs. 
  # Then create a new dataframe with the colnames of viraltop100 and populate with the raw data info for remaining raw data species and 0/NA for sections not present.
  # finally rbind the finished table. 
  
  # Then we figure out confidence metrics
  
  
  
  for (i in c(1:nrow(Viraltop100))) {
    
    
    grep(pattern=Viraltop100$species[i],x=viral_raws_complexity_summary$species) -> idx
    
    if (length(idx) > 0) {
      
      Viraltop100$additional_raw_reads_assigned[i] <- "yes"
      Viraltop100$number_of_raw_reads_assigned[i] <- viral_raws_complexity_summary$count[idx]
      Viraltop100$mean_complexity_of_raw_reads[i] <- viral_raws_complexity_summary$mean_complexity[idx]
      Viraltop100$mean_identity_of_raw_reads[i] <- viral_raws_complexity_summary$mean_identity[idx]
      Viraltop100$mean_aligned_length_of_raw_reads[i] <- viral_raws_complexity_summary$mean_aligned_length[idx]
      
      
    }
    
  }
  
}



viral_raws_complexity_summary <- subset(viral_raws_complexity_summary,viral_raws_complexity_summary$superkingdom=="Viruses")

viral_raws_complexity_summary$subspecies <- as.character(viral_raws_complexity_summary$subspecies)
Viraltop100$subspecies <- as.character(Viraltop100$subspecies)

# First check if the raws is empty.
# second check if the contigs is empty. (heirarchically)

if (file.info(viralcomplexity)$size >= 350) {
  remaining_raws <- anti_join(viral_raws_complexity_summary, Viraltop100, by = "subspecies")
  
  if(nrow(Combined_assigned_contigs_viruses_only) >=1) {
    
    viral_raws_reads_compatability <- as.data.frame(matrix(ncol=ncol(Viraltop100), nrow=nrow(remaining_raws)))
    colnames(viral_raws_reads_compatability) <- colnames(Viraltop100)
    BOTHMISSING <- "NO"
  }
  if(nrow(Combined_assigned_contigs_viruses_only) <1) {
    
    viral_raws_reads_compatability <- as.data.frame(matrix(ncol=21, nrow=nrow(remaining_raws)))
    colnames(viral_raws_reads_compatability) <- c("reads_assigned_through_contigs","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident","max_percent_ident","length","number_contigs_assigned","mean_complexity_of_contigs","additional_raw_reads_assigned","number_of_raw_reads_assigned","mean_complexity_of_raw_reads","mean_identity_of_raw_reads","mean_aligned_length_of_raw_reads")
    BOTHMISSING <- "NO"
  }
  
  
}

if (file.info(viralcomplexity)$size < 350) {
  
  if(nrow(Combined_assigned_contigs_viruses_only) >=1) {
    
    viral_raws_reads_compatability <- as.data.frame(matrix(ncol=ncol(Viraltop100), nrow=0))
    colnames(viral_raws_reads_compatability) <- colnames(Viraltop100)
    BOTHMISSING <- "NO"
  }
  if(nrow(Combined_assigned_contigs_viruses_only) <1) {
    
    viral_raws_reads_compatability <- as.data.frame(matrix(ncol=21, nrow=0))
    colnames(viral_raws_reads_compatability) <- c("reads_assigned_through_contigs","superkingdom","phylum","class","order","family","genus","species","subspecies","average_percent_ident","min_percent_ident","max_percent_ident","max_percent_ident","length","number_contigs_assigned","mean_complexity_of_contigs","additional_raw_reads_assigned","number_of_raw_reads_assigned","mean_complexity_of_raw_reads","mean_identity_of_raw_reads","mean_aligned_length_of_raw_reads")
    BOTHMISSING <- "YES"
  }
  
}


if (BOTHMISSING =="NO") {
  
  
  if (file.info(viralcomplexity)$size >= 350) {
    
    viral_raws_reads_compatability$reads_assigned_through_contigs <- 0
    viral_raws_reads_compatability$superkingdom <- remaining_raws$superkingdom
    viral_raws_reads_compatability$phylum <- remaining_raws$phylum
    viral_raws_reads_compatability$class <- remaining_raws$class
    viral_raws_reads_compatability$order <- remaining_raws$order
    viral_raws_reads_compatability$family <- remaining_raws$family
    viral_raws_reads_compatability$contigs_assigned <- "no"
    viral_raws_reads_compatability$genus <- remaining_raws$genus
    viral_raws_reads_compatability$species <- remaining_raws$species
    viral_raws_reads_compatability$subspecies <- remaining_raws$subspecies
    viral_raws_reads_compatability$average_percent_ident <- 0
    viral_raws_reads_compatability$min_percent_ident <- 0
    viral_raws_reads_compatability$max_percent_ident <- 0
    viral_raws_reads_compatability$length <- 0
    viral_raws_reads_compatability$number_contigs_assigned <- 0
    viral_raws_reads_compatability$mean_complexity_of_contigs <- 0
    viral_raws_reads_compatability$additional_raw_reads_assigned <- "yes" 
    viral_raws_reads_compatability$number_of_raw_reads_assigned <- remaining_raws$count
    viral_raws_reads_compatability$mean_complexity_of_raw_reads <- remaining_raws$mean_complexity
    viral_raws_reads_compatability$mean_identity_of_raw_reads <- remaining_raws$mean_identity
    viral_raws_reads_compatability$mean_aligned_length_of_raw_reads <- remaining_raws$mean_aligned_length
  }
  
  
  if (file.info(viralcomplexity)$size < 350) {
    
    print("no raw reads assigned to viruses")
  }
  

  
  
  if(nrow(Combined_assigned_contigs_viruses_only) >=1) {
    
    
    if (file.info(viralcomplexity)$size >= 350) {
      
      virus_all <- rbind(Viraltop100,viral_raws_reads_compatability)
      
    }
    
    
    
    if (file.info(viralcomplexity)$size < 350) {
      
      virus_all <- Viraltop100
      
    }
    
    
    
  }
  
  if(nrow(Combined_assigned_contigs_viruses_only) <1) {
    
    virus_all <- viral_raws_reads_compatability
    
  }
  
  
  virus_all$total_reads_assigned <- virus_all$reads_assigned_through_contigs + virus_all$number_of_raw_reads_assigned
  
  virus_all$classification_level_estimate <- ""
  virus_all$confidence_of_assignment_as_virus <- ""
  
  
  
  
  virus_all$reads_assigned_through_contigs <- as.numeric(virus_all$reads_assigned_through_contigs)
  virus_all$number_of_raw_reads_assigned <- as.numeric(virus_all$number_of_raw_reads_assigned)
  virus_all$mean_complexity_of_contigs <- as.numeric(virus_all$mean_complexity_of_contigs)
  virus_all$mean_complexity_of_raw_reads <- as.numeric(virus_all$mean_complexity_of_raw_reads)
  virus_all$average_percent_ident <- as.numeric(virus_all$average_percent_ident)
  virus_all$mean_identity_of_raw_reads <- as.numeric(virus_all$mean_identity_of_raw_reads)
  virus_all$length <- as.numeric(virus_all$length)
  

  # Rename 'staxidreduced' to 'taxid' in Combined_assigned_contigs_viruses_only so that left joins work 
  Combined_assigned_contigs_viruses_only <- Combined_assigned_contigs_viruses_only %>%
    rename(taxid = staxidreduced) %>%
    mutate(taxid = as.character(taxid))  # Convert taxid to character

  # Ensure taxid column exists in virus_all and is character
  colnames(virus_all) <- make.unique(colnames(virus_all))

  virus_all <- virus_all %>%
    mutate(taxid = as.character(NA))

  # Make sure everything is the same object class
  remaining_raws <- remaining_raws %>%
    mutate(taxid = as.character(taxid))

  # left joins grab the taxids from both the contig and the raws dfs
  virus_all <- virus_all %>%
    left_join(Combined_assigned_contigs_viruses_only %>% select(subspecies, taxid), by = "subspecies", suffix = c("", "_contigs")) %>%
    left_join(remaining_raws %>% select(subspecies, taxid), by = "subspecies", suffix = c("", "_raws")) %>%
    mutate(taxid = coalesce(taxid, taxid_contigs, taxid_raws)) %>%  # Fill taxid in priority order
    select(-taxid_contigs, -taxid_raws)  # Remove extra columns

  virus_all <- virus_all %>%
    distinct(subspecies, .keep_all = TRUE)

  virus_all <- virus_all %>%
    mutate(taxid = sub("[:;].*", "", taxid))


  virus_all <- virus_all %>%
    arrange(desc(total_reads_assigned))





  colourise <- function(value, column) {
    if (column == "reads_assigned_through_contigs") {
      if (value >= 50) {
        return("green")
      } else if (value >= 25) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else {
        return("red")
      }
    } else if (column == "total_reads_assigned") {
      if (value >= 500) {
        return("green")
      } else if (value >= 1 && value <=20) {
        return("red")
      } else {
        return(NA)
      }
    } else if (column == "number_of_raw_reads_assigned") {
      if (value >= 50) {
        return("green")
      } else if (value >= 25) {
        return("orange")
      }else if (value ==0) {
        return(NA)
      } else {
        return("red")
      }
    } else if (column == "mean_complexity_of_contigs") {
      if (value >= 45) {
        return("green")
      } else if (value >= 33) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else {
        return("red")
      }
    } else if (column == "mean_complexity_of_raw_reads") {
      if (value >= 40) {
        return("green")
      } else if (value >= 33) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else {
        return("red")
      }
    } else if (column == "average_percent_ident") {
      if (value >= 90) {
        return("green")
      } else if (value >= 80) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else {
        return("blue")
      }
    } else if (column == "mean_identity_of_raw_reads") {
      if (value >= 90) {
        return("green")
      } else if (value >= 80) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else if (value <= 80){
        return("blue")
      }
    } else if (column == "length") {
      if (value >= 600) {
        return("green")
      } else if (value >= 300) {
        return("orange")
      } else if (value ==0) {
        return(NA)
      } else {
        return("red")
      }
    }    
  }
  
  coloured_data <- sapply(names(virus_all), function(col) {
    sapply(virus_all[[col]], function(x) {
      paste0("<span style='color:", colourise(x, col), "'>", x, "</span>")
    })
  }, simplify = FALSE)
  
  
  coloured_data <- bind_cols(coloured_data)
  
  #coloured_data <- as.data.frame(coloured_data)
  
  
  for (i in c(1:nrow(coloured_data))) {
    
    greens <- sum(grepl(pattern="color:green",x=coloured_data[i,c(1, 13, 15, 18, 19, 22)]))
    oranges <- sum(grepl(pattern="color:orange",x=coloured_data[i,c(1, 13, 15, 18, 19, 22)]))
    reds <- sum(grepl(pattern="color:red",x=coloured_data[i,c(1, 13, 15, 18, 19, 20)]))
    
    
    
    if(greens==1 && reds<=1) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Medium"
      
      
    }
    
    
    
    if(greens>=2 &&reds<=1) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "High"
      
      
    }
    
    
    if(greens==2 &&reds==1) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Medium"
      
      
    }
    
    if(greens>=2 &&reds>=2) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Medium"
      
      
    }
    
    if(greens<2 &&reds>=2) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Low"
      
      
    }
    
    if(greens==0 && oranges>=2) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Low"
      
      
    }
    
    
    if(greens==0 && oranges>=1 && reds >=1) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Low"
      
      
    }
    
    if(greens==2 && reds == 0 && oranges == 0) {
      
      coloured_data$confidence_of_assignment_as_virus[i] <- "Medium-High"
      
      
    }
    
    
    
    greens2 <- sum(grepl(pattern="color:green",x=coloured_data[i,c(10, 13, 20)]))
    oranges2 <- sum(grepl(pattern="color:orange",x=coloured_data[i,c(10, 13, 20)]))
    reds2 <- sum(grepl(pattern="color:red",x=coloured_data[i,c(10, 13, 20)]))  
    blues2 <- sum(grepl(pattern="color:blue",x=coloured_data[i,c(10, 13, 20)]))  
    totcols <- sum(greens2 + oranges2 + reds2 + blues2)
    if(greens2 >=2) {
      
      coloured_data$classification_level_estimate <- "Likely same species"
    }
    
    if(greens2 >=2 && oranges2 ==0 && reds2 ==0 && blues2 == 0) {
      
      coloured_data$classification_level_estimate[i] <- "Likely same species"
    }
    
    if(greens2 ==1 && oranges2 ==0 && reds2 ==0 && blues2 == 0) {
      
      coloured_data$classification_level_estimate[i] <- "Likely same species but recommend double checking"
    }
    
    if(greens2 >=2 && oranges2 ==1 && reds2 ==0 && blues2 == 0) {
      
      coloured_data$classification_level_estimate[i] <- "Likely same species"
    }
    
    if(greens2 >=1 && oranges2 >=2 ) {
      
      coloured_data$classification_level_estimate[i] <- "Possible closely related species"
    }
    
    
    if(reds2 ==1) {
      
      coloured_data$classification_level_estimate[i] <- "Possible closely related species"
    }
    
    if(reds2 >=2) {
      
      coloured_data$classification_level_estimate[i] <- "Possible closely related species"
    }
    
    if(blues2 >=1) {
      
      coloured_data$classification_level_estimate[i] <- "Likely different species"
    }
    
    if(totcols <=1 && greens2 ==0) {
      
      coloured_data$classification_level_estimate[i] <- "Insufficient data to reliably infer"
    }


  }
  

 
  
  virus_all$classification_level_estimate <- coloured_data$classification_level_estimate
  virus_all$confidence_of_assignment_as_virus <- coloured_data$confidence_of_assignment_as_virus
  
  confidence_order <- c("High", "Medium-High", "Medium", "Low")
  



  virus_all$confidence_of_assignment_as_virus <- factor(virus_all$confidence_of_assignment_as_virus, levels = confidence_order)

  




  virus_all <- virus_all %>%
    select(
      reads_assigned_through_contigs, superkingdom, phylum, class, order, family, genus, species, subspecies, 
      taxid, everything()
    )


  # For diverged read and contig detection I want to make sure a good numer of reads/contigs were detected already (set as 10) to minimise FP, plus I want to subset to only 
  # species where average idents are below 95% as working with something about that suggests that the detected species is a species/strain match to the reference and therefore
  # introducing divergent reads analysis is disproportionately likely to generate false positives.   
  virus_all2 <- subset(virus_all,virus_all$total_reads_assigned>10 & (virus_all$average_percent_ident <95 | virus_all$mean_identity_of_raw_reads <95))

  
  coloured_data <- coloured_data %>%
    select(
      reads_assigned_through_contigs, superkingdom, phylum, class, order, family, genus, species, subspecies, 
      taxid, everything()
    )


  
  coloured_datatable <- datatable(coloured_data, escape = FALSE)
  
  coloured_datatable <- datatable(coloured_data, escape = FALSE, callback = JS("return table;"))
  
  
  
  # Completed generating full viral identification table. (But may come back and incorporate a viral trees component)
  
  # Now to extract all reads and all contigs assosciated with each species by building species lists. 
  
  #start with Raw reads. Don't forget to also have the raws assigned to each of the contigs from step 9.
  # Three dataframes needed here. 
  #1. the virus_all for species names (generated here)
  #2. Diamondrawhitsvirus for read names of raw assigned virus reads (file from R env imported as Diamondrawshits then filtered down to just virus)
  #3. raw reads assigned to each contig (will be in file from R env, added as of 5/4/24) # confirmed in new pipeline and outputting correctly 
  #4. allassignedfreqsvirus for the list of contigs assigned to each viral species.  # confirmed in new pipeline and outputting correctly 
  
  #1. virus_all
  #2. Diamondrawshitsvirus
  #3. reads_contigs_virus
  #4. allassignedfreqsvirus
  #5. the full raws assignment table to pick up extra kraken specific reads

rawsassigntable$Sample <- NAMES
colnames(rawsassigntable) <- colnames(Diamondrawshits)

rawsassigntablevirus <- subset(rawsassigntable,rawsassigntable$superkingdom=="Viruses")  


  Diamondrawhitsvirus <- subset(Diamondrawshits, Diamondrawshits$superkingdom=="Viruses")
  Diamondrawhitsvirusorig <- Diamondrawhitsvirus

Diamondrawhitsvirus2 <- rbind(rawsassigntablevirus,Diamondrawhitsvirus)
Diamondrawhitsvirus <- Diamondrawhitsvirus2

  combined_virus_reads <- as.data.frame(matrix(nrow=0,ncol=2))
  colnames(combined_virus_reads) <- c("Reads","Species")
  combined_virus_contigs <- as.data.frame(matrix(nrow=0,ncol=2))
  colnames(combined_virus_contigs) <- c("Contigs","Species")
  virus_all_uniqsp <- virus_all %>% distinct(species, .keep_all = TRUE)


for (i in c(1:nrow(virus_all_uniqsp))) {
  
  
  # Start with contigs
  workingcontigs <- subset(allassignedfreqsvirus,allassignedfreqsvirus$species==virus_all_uniqsp$species[i])
  
  contigsdf <- as.data.frame(matrix(nrow=nrow(workingcontigs),ncol=2))
  contigsdf$V1 <- workingcontigs[,1]
  
  if (nrow(contigsdf) >=1){
    contigsdf$V2 <- virus_all_uniqsp$species[i]
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
  reads_contigs <- as.data.frame(matrix(nrow=length(idx2),ncol=2))
  reads_contigs$V1 <- reads_contigs_virus[idx2,1]
  
  if (nrow(reads_contigs) >=1){
    reads_contigs$V2 <- virus_all_uniqsp$species[i]
  }
  
  colnames(reads_contigs) <- c("Reads","Species")
  
  # Now do raws
  grep(pattern=virus_all_uniqsp$species[i],x=Diamondrawhitsvirus$species) ->idx
  
  reads_raws <- as.data.frame(matrix(nrow=length(idx),ncol=2))
  reads_raws$V1 <- Diamondrawhitsvirus[idx,1]
  
  if (nrow(reads_raws) >=1){
    reads_raws$V2 <- virus_all_uniqsp$species[i]
  }
  colnames(reads_raws) <- c("Reads","Species")
  
  reads_combined_singlevirus <- rbind(reads_contigs,reads_raws)
  
  combined_virus_reads <- rbind(combined_virus_reads,reads_combined_singlevirus)
  
  combined_virus_contigs <- rbind(combined_virus_contigs,contigsdf)
  
}  

  # So there is a minor issue with the above (more about pair read data in general). This returns duplicates of each read in the pair
  #as the read 2:N:0 or 1:N:0 is cut off
  # This is solved by removing duplicates but in the scenario that one read of a pair assigned as virus but not the other, this will carry over the non virus read to
  # But the question is, is that bad? The pairs should be pairs and part of the same sequence. Further, if I took only the single reads then the fq files wouldn't be
  # paired anymore and couldn't be analysed further manually as order would be disrupted. 
  # I didn't realise what the problem was earlier but the reason the raws still have the cassava formatting but the finals don't is because
  # sam formats fail to keep formatted headers with a whitespace in them. After playing around I have managed to reattach the reads for steps 3-5 
  # but there may be other instances of this issue what I haven't yet picked up on in the pipeline. (I will keep an eye out but as above it shouldn't actually
  # matter, plus, as described above, returning uneven R1s and R2s messes with a lot of programs.)
  # I may need to incorporate an unaligned reads.fastq file?  
  
  
  length(unique(combined_virus_reads$Reads))
  combined_virus_reads_unique <- unique(combined_virus_reads)
  # Note all contigs were unique already though
  combined_virus_contigs_unique <- unique(combined_virus_contigs)
  
  combined_virus_reads_unique$Species <- gsub(" ", "_", combined_virus_reads_unique$Species)
  combined_virus_contigs_unique$Species <- gsub(" ", "_", combined_virus_contigs_unique$Species)
  
  # Final step. I need to iterate a save loop which will save a file with output/NAME/species_reads_title
  # and output/NAME/species_contigs_title
  
  
  # Define the directory to save the text files
  # The input params for output directory includes the NAMES already 
  directory <- resultspathcompile
  
  #directory <- paste0(getwd(),"/",NAMES,"/")
  
  # Create the directory if it doesn't exist
  if (!dir.exists(directory)) {
    dir.create(directory,recursive = TRUE)
  }
  
  # Split the dataframe into a list of dataframes based on the values in the category column
  split_dfraws <- split(combined_virus_reads_unique, combined_virus_reads_unique$Species)

# Some names have slashes in them which mess with the directory creation so change them to underscores
  names(split_dfraws) <- gsub("/", "_", names(split_dfraws))

  # Save each dataframe in the list as a separate text file
  for (category_value in names(split_dfraws)) {
    file_name <- paste0(directory, category_value, "/",category_value,"_read_names.txt")
    Readnames <- split_dfraws[[category_value]][, 1]
    
    dir.create(paste0(directory, category_value, "/"))
    
    write.table(Readnames, file = file_name, sep = "\t", row.names = FALSE,col.names=FALSE,quote = FALSE)
  }
  

  # repeat for contigs
  split_dfcontigs <- split(combined_virus_contigs_unique, combined_virus_contigs_unique$Species)
names(split_dfcontigs) <- gsub("/", "_", names(split_dfcontigs))
  
  # Save each dataframe in the list as a separate text file
  for (category_value in names(split_dfcontigs)) {
    file_name <- paste0(directory, category_value, "/",category_value,"_contig_names.txt")
    contignames <- split_dfcontigs[[category_value]][, 1]
    write.table(split_dfcontigs[[category_value]], file = file_name, sep = "\t", row.names = FALSE,col.names=FALSE,quote = FALSE)
  }
  
  
  
  
  write.csv(virus_all, file=paste0(outtablespathcompile,NAMES,"_virusall_sums.csv"), row.names = FALSE)

  write.table(virus_all2$taxid, file = paste0(outtablespathcompile, NAMES, "_detected_viral_taxids.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
  
  datatable(coloured_data, escape = FALSE) %>%
    saveWidget(file = paste0(outtablespathcompile,NAMES,"_virusall_datatable.html"), selfcontained = TRUE)
  
}


if (BOTHMISSING =="YES") {
  
  file.create(file = paste0(outtablespathcompile,NAMES,"_virusall_datatable.html"))
  file.create(file=paste0(outtablespathcompile,NAMES,"_virusall_sums.csv"))
  file.create(file=paste0(outtablespathcompile,NAMES,"_detected_viral_taxids.csv"))
  
}

# So need to update the arg pass values so they point correctly, then create the rule in the analyse_raws rule smk
# Then need to create another rule to use zcat | grep to extract all reads and contigs associated with each species of virus. 
# For the reads, grep -A 4 will work as the reads are all 4line. For the contigs I hope it will work, but if not I may need to use something
#like seqkt to do it

#save.image(file=paste0(outtablespathcompile,NAMES,"_output_results_R_env.Rdata"))

