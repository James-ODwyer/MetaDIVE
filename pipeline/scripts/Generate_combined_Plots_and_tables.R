# Script to compile multiple sample summary results into 1, combined tables, 2 create summary barplots which show where everything is assigned for all reads
# and 3, for specific species

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
library(viridis)
library(hrbrthemes)
library(reshape2)
library(scales)



if (!require("d3treeR", character.only = TRUE)) {
  devtools::install_github("timelyportfolio/d3treeR",upgrade="never",dependencies=FALSE)
  library("d3treeR")
} else {
  library("d3treeR")
}


parser <- ArgumentParser(description= 'Summarising results filtering and assembly')

parser$add_argument('--inputpath', '-a', help= 'input directory for summary results to compile')
parser$add_argument('--outputpath', '-b', help= 'Output directory to save figures and combined tables')
parser$add_argument('--programdir', '-c', help= 'Base program directory')
parser$add_argument('--dohostdetect', '-d', help= 'whether detect and remove host Eukaryotic species was undertaken')
parser$add_argument('--dodiamondraws', '-e', help= 'Whether raw reads were blasted in Diamondx')

xargs<- parser$parse_args()

# define base parameters and parameter variables for which analyses were run
basepath <- xargs$programdir
resultspath <- xargs$outputpath
inputdir <- xargs$inputpath
outtablespath <- paste0(basepath,resultspath)
intablespath <- paste0(basepath,inputdir)
dohostdetect <- xargs$dohostdetect
dodiamondraws <- xargs$dodiamondraws


samplereadssummaryfiles <- list.files(path = intablespath, pattern = "assignment_reads_for_plot_generation", all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE)


Allreadssummaryresultstable <- list()

for (i in c(1:length(samplereadssummaryfiles))) {

  Resultstable <- read.table(samplereadssummaryfiles[i],header = TRUE,sep = "\t")
  sample<- Resultstable$Sample
  
  Resultstable <- Resultstable[,-11]
  
  Resultstablet <- t(Resultstable)
  
  Resultstablet[,1]
  
  Resultstableformatted <- as.data.frame(matrix(nrow=10,ncol=3))
  
  Resultstableformatted$V1 <- rownames(Resultstablet)
  Resultstableformatted$V2 <- Resultstablet[,1]
  Resultstableformatted$V3 <- sample
  
  colnames(Resultstableformatted) <- c("Filtering_step","Reads","Sample")
  
  Allreadssummaryresultstable[[i]] <- Resultstableformatted
}

Resultssummaryreadsassignments <-bind_rows(Allreadssummaryresultstable)


#Resultssummaryreadsassignments$Filtering_step  <- factor(Resultssummaryreadsassignments$Filtering_step,levels = c("Low.quality.reads", "PhiX.contamination", "CO1.hits", "LSU.hits", "SSU.hits", "Host.species", "Eukaryotes", "Bacteria", "Viruses", "Unassigned.reads"))
Resultssummaryreadsassignments$Filtering_step  <- factor(Resultssummaryreadsassignments$Filtering_step,levels = c("Unassigned.reads", "Viruses", "Bacteria", "Eukaryotes", "Host.species", "SSU.hits", "LSU.hits", "CO1.hits", "PhiX.contamination", "Low.quality.reads"))

Resultssummaryreadsassignments$Reads[is.na(Resultssummaryreadsassignments$Reads)] <- 0

plot <- ggplot(Resultssummaryreadsassignments, aes(x=Sample, y=Reads, fill=Filtering_step)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Breakdown of total reads assignments per sample") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7),axis.text.x = element_text(size=4.5)) +
  theme(legend.title = element_text(size=6,), legend.text = element_text(size=4.5)) +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=4.5)) +
  theme(plot.title =  element_text(size=9,hjust = 0.4)) +
  scale_y_continuous(labels = scales::number_format()) +
  coord_flip() +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot 
ggsave(paste0(outtablespath,"summary_all_reads_assigned_and_filtered.pdf"), dpi = 300)

plot 
ggsave(paste0(outtablespath,"summary_all_reads_assigned_and_filtered.png"), dpi = 300)


# Generate large combined reads filtering summary table
summary_reads_filteringtables <- list.files(path = intablespath, pattern = "summary_reads_filtering", all.files = FALSE,
                                      full.names = TRUE, recursive = FALSE,
                                      ignore.case = FALSE, include.dirs = FALSE)


combined_summary_reads_filteringtables <- as.data.frame(matrix(nrow=length(summary_reads_filteringtables),ncol=13))

for (i in c(1:length(summary_reads_filteringtables))) {
  
  Resultstable <- read.table(summary_reads_filteringtables[i],header = TRUE,sep = "\t")

  
  combined_summary_reads_filteringtables[i,] <- Resultstable
  
  
}

colnames(combined_summary_reads_filteringtables) <- colnames(Resultstable)


write.table(combined_summary_reads_filteringtables,file=(paste0(outtablespath,"Combined_summary_reads_filtering.txt")),sep="\t",row.names=FALSE)


# Generate large combined contigs summary table
summary_contigs_filteringtables <- list.files(path = intablespath, pattern = "summarycontigs_assembly_values", all.files = FALSE,
                                            full.names = TRUE, recursive = FALSE,
                                            ignore.case = FALSE, include.dirs = FALSE)


combined_summary_contigs_filteringtables <- as.data.frame(matrix(nrow=length(summary_contigs_filteringtables),ncol=17))

for (i in c(1:length(summary_contigs_filteringtables))) {
  
  Resultstable <- read.table(summary_contigs_filteringtables[i],header = TRUE,sep = "\t")
  
  
  combined_summary_contigs_filteringtables[i,] <- Resultstable
  
  
}

colnames(combined_summary_contigs_filteringtables) <- colnames(Resultstable)

write.table(combined_summary_contigs_filteringtables,file=(paste0(outtablespath,"Combined_summary_contigs_assiging.txt")),sep="\t",row.names=FALSE)

# Now bring in the Contig assigned read files to generate figures for top returned species. 

# Binary yes no to whether to read in the no host freq summary or the full freq summary

# It can read on empty files, but then there is no sample name which was stored in the files otherwise

# Top 5 of each superkingdom as separate lists per sample. run files to create and assign each. Then rbind list values into combined data frame for 
# combined figure generation.


if (dohostdetect=="yes") {
  
 
   summaryreturnedhitsfromcontigs <- list.files(path = intablespath, pattern = "summarycontighits_assigned_no_host_assembly", all.files = FALSE,
                                               full.names = TRUE, recursive = FALSE,
                                               ignore.case = FALSE, include.dirs = FALSE)
  
  
  
  
}

if (dohostdetect=="no") {
  
  
  summaryreturnedhitsfromcontigs <- list.files(path = intablespath, pattern = "summarycontighits_assigned_assembly", all.files = FALSE,
                                               full.names = TRUE, recursive = FALSE,
                                               ignore.case = FALSE, include.dirs = FALSE)
  
}


Euklist <- list()
bactlist <- list()
virlist <- list()



for (i in c(1:length(summaryreturnedhitsfromcontigs))) {
  
  Resultstable <- read.table(summaryreturnedhitsfromcontigs[i],header = TRUE,sep = "\t")
  sampname <- gsub(".*99_SUMMARY_RESULTS/",x=summaryreturnedhitsfromcontigs[i],replacement="")
  sampname <- gsub("_summarycontighits.*",x=sampname,replacement="")
  sampname <- gsub("/",x=sampname,replacement="")
  Resultstable <- subset(Resultstable, !(is.na(Resultstable$species)))
  # Need to remove all NA rows from data
  
  ResultstableEuk <- subset(Resultstable,Resultstable$superkingdom=="Eukaryota")
  
  uniquespecies <- unique(ResultstableEuk$species)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=6))
  colnames(speciescounts) <- c("Species", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(ResultstableEuk)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(ResultstableEuk,ResultstableEuk$species==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- sum(resultstablesubset$freq)
    speciescounts[j,3] <- mean(resultstablesubset$percentident)
    speciescounts[j,4] <- max(resultstablesubset$percentident)
    speciescounts[j,5] <- min(resultstablesubset$percentident)
    speciescounts[j,6] <- resultstablesubset$Sample[1]
    
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }
  Euklist[[i]] <- speciescounts
  
}
  
  if (nrow(ResultstableEuk)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    Euklist[[i]] <- speciescounts
  }
  
  
  
  
  Resultstablebacteria <- subset(Resultstable,Resultstable$superkingdom=="Bacteria")
  
  
  uniquespecies <- unique(Resultstablebacteria$species)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=6))
  colnames(speciescounts) <- c("Species", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(Resultstablebacteria)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(Resultstablebacteria,Resultstablebacteria$species==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- sum(resultstablesubset$freq)
    speciescounts[j,3] <- mean(resultstablesubset$percentident)
    speciescounts[j,4] <- max(resultstablesubset$percentident)
    speciescounts[j,5] <- min(resultstablesubset$percentident)
    speciescounts[j,6] <- resultstablesubset$Sample[1]
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }
  bactlist[[i]] <- speciescounts
  
}
  
  if (nrow(Resultstablebacteria)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    bactlist[[i]] <- speciescounts
  }
  
  ResultstableVirus <- subset(Resultstable,Resultstable$superkingdom=="Viruses")
  
  uniquespecies <- unique(ResultstableVirus$species)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=6))
  colnames(speciescounts) <- c("Species", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(ResultstableVirus)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(ResultstableVirus,ResultstableVirus$species==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- sum(resultstablesubset$freq)
    speciescounts[j,3] <- mean(resultstablesubset$percentident)
    speciescounts[j,4] <- max(resultstablesubset$percentident)
    speciescounts[j,5] <- min(resultstablesubset$percentident)
    speciescounts[j,6] <- resultstablesubset$Sample[1]
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }
  virlist[[i]] <- speciescounts
  }
  
  if (nrow(ResultstableVirus)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    virlist[[i]] <- speciescounts
  }
  
  
  
}
# Combine lists
Viruses <- bind_rows(virlist)
Bacteria <- bind_rows(bactlist)
Eukaryotes <- bind_rows(Euklist)

# Now create top 10s based on aggregate reads assigned.

# Viruses
Viruses_sums <- Viruses[order(-Viruses$`Reads assigned`),]

#Viruses_sums <- aggregate(Viruses,by = Viruses$`Reads assigned`,FUN = sum())
Viruses_sums <- aggregate((Viruses$`Reads assigned`),by=list(Viruses$Species),sum)
Viruses_sums <- Viruses_sums[order(-Viruses_sums$`x`),]

if (nrow(Viruses_sums)>=10) {
  Viruses_sumstop10 <- Viruses_sums[1:10,]
}
if (nrow(Viruses_sums)<10) {
  Viruses_sumstop10 <- Viruses_sums
}

Viruses_top10 <- Viruses %>% 
  filter(Species %in% Viruses_sumstop10$Group.1)

Viruses_sumstop10$Group.1


# Bacteria
Bacteria_sums <- Bacteria[order(-Bacteria$`Reads assigned`),]

#Viruses_sums <- aggregate(Viruses,by = Viruses$`Reads assigned`,FUN = sum())
Bacteria_sums <- aggregate((Bacteria$`Reads assigned`),by=list(Bacteria$Species),sum)
Bacteria_sums <- Bacteria_sums[order(-Bacteria_sums$`x`),]

if (nrow(Bacteria_sums)>=10) {
  Bacteria_sumstop10 <- Bacteria_sums[1:10,]
}
if (nrow(Bacteria_sums)<10) {
  Bacteria_sumstop10 <- Bacteria_sums
}

Bacteria_top10 <- Bacteria %>% 
  filter(Species %in% Bacteria_sumstop10$Group.1)

Bacteria_top10$Group.1


# Eukaryotes
Eukaryotes_sums <- Eukaryotes[order(-Eukaryotes$`Reads assigned`),]

#Viruses_sums <- aggregate(Viruses,by = Viruses$`Reads assigned`,FUN = sum())
Eukaryotes_sums <- aggregate((Eukaryotes$`Reads assigned`),by=list(Eukaryotes$Species),sum)
Eukaryotes_sums <- Eukaryotes_sums[order(-Eukaryotes_sums$`x`),]

if (nrow(Eukaryotes_sums)>=10) {
  Eukaryotes_sumstop10 <- Eukaryotes_sums[1:10,]
}
if (nrow(Eukaryotes_sums)<10) {
  Eukaryotes_sumstop10 <- Eukaryotes_sums
}

Eukaryotes_top10 <- Eukaryotes %>% 
  filter(Species %in% Eukaryotes_sumstop10$Group.1)

Eukaryotes_top10$Group.1






# ggplot graphs 

Viruses_top10$`Reads assigned` <- as.numeric(as.character(Viruses_top10$`Reads assigned`))
Bacteria_top10$`Reads assigned` <- as.numeric(as.character(Bacteria_top10$`Reads assigned`))
Eukaryotes_top10$`Reads assigned` <- as.numeric(as.character(Eukaryotes_top10$`Reads assigned`))


# ggplot graphs 


plot <- ggplot(Viruses_top10, aes(x=Sample, y=`Reads assigned`, fill=Species)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top 10 virus species reads/contigs \n  were assigned to") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 18)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 16),axis.text.x = element_text(size=14)) +
  theme(legend.title = element_text(size=12), legend.text = element_text(size=4.0)) +
  theme(legend.key.size = unit(0.32,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=16.5)) +
  theme(plot.title =  element_text(size=19,hjust = 0.4)) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(labels = scales::number_format()) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot

ggsave(paste0(outtablespath,"Summary_barplot_top_virus_species_from_contigs.pdf"), dpi = 300)

plot 

ggsave(paste0(outtablespath,"Summary_barplot_top_virus_species_from_contigs.png"), dpi = 300)

plot <- ggplot(Bacteria_top10, aes(x=Sample, y=`Reads assigned`, fill=Species)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top 10 Bacterial species reads/contigs \n  were assigned to") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 18)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 16),axis.text.x = element_text(size=14)) +
  theme(legend.title = element_text(size=12), legend.text = element_text(size=4.0)) +
  theme(legend.key.size = unit(0.33,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=16.5)) +
  theme(plot.title =  element_text(size=20,hjust = 0.4)) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(labels = scales::number_format()) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot


ggsave(paste0(outtablespath,"summarybarplot_contigs_bacterial_species.pdf"), dpi = 300)

plot 

ggsave(paste0(outtablespath,"summarybarplot_contigs_bacterial_species.png"), dpi = 300)


plot <- ggplot(Eukaryotes_top10, aes(x=Sample, y=`Reads assigned`, fill=Species)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top 10 Eukaryote species reads/contigs \n  were assigned to") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 18)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 16),axis.text.x = element_text(size=14)) +
  theme(legend.title = element_text(size=12), legend.text = element_text(size=4.0)) +
  theme(legend.key.size = unit(0.33,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=16.5)) +
  theme(plot.title =  element_text(size=20,hjust = 0.4)) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(labels = scales::number_format()) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot


ggsave(paste0(outtablespath,"summarybarplot_contigs_Eukaryote_species.pdf"), dpi = 300)

plot 

ggsave(paste0(outtablespath,"summarybarplot_contigs_Eukaryote_species.png"), dpi = 300)


# Now to repeat the above but for raw Diamond hits and at the family level to compensate for the lower accuracy likely found.

# Make sure it is only if do diamond raws 

if (dodiamondraws=="yes") {


summary_diamondrawstables <- list.files(path = intablespath, pattern = "_raw_diamond_hits.txt", all.files = FALSE,
                                              full.names = TRUE, recursive = FALSE,
                                              ignore.case = FALSE, include.dirs = FALSE)


combined_summary_diamondrawstables <- as.data.frame(matrix(nrow=length(summary_diamondrawstables),ncol=17))

# So Haven't switched these over to top 10 yet?


Euklistdiamondraws <- list()
bactlistdiamondraws  <- list()
virlistdiamondraws  <- list()


# After doing some prelim testing on most effective cut offs for pw ident in diamond raws. 
# Even at 75% pw identity min, 6 different carnivor families other than canid was found. This suggests that 
# 75% still is overly vague (for canids at least) and so anything lower would likely drive more off assignments. I don't want to 
# go much higher as almost half of reads are already filtered and it is a general just observe what is likely in there not actually classify anything

for (i in c(1:length(summary_diamondrawstables))) {
  
  Resultstable <- read.table(summary_diamondrawstables[i],header = TRUE,sep = "\t")
  sampname <- gsub(".*99/",x=summary_diamondrawstables[i],replacement="")
  sampname <- gsub("_raw_diamond_hits.*",x=sampname,replacement="")
  Resultstable <- subset(Resultstable, !(is.na(Resultstable$family)))
  # Need to remove all NA rows from data
  Resultstable <- subset(Resultstable,Resultstable$pident>=75)
  ResultstableEuk <- subset(Resultstable,Resultstable$superkingdom=="Eukaryota")
  
  uniquefamilies <- unique(ResultstableEuk$family)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquefamilies), ncol=6))
  colnames(speciescounts) <- c("Family", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(ResultstableEuk)>=1) { 
    for ( j in c(1:length(uniquefamilies))) {
      
      resultstablesubset <- subset(ResultstableEuk,ResultstableEuk$family==uniquefamilies[j])
      
      speciescounts[j,1] <- resultstablesubset$family[1]
      speciescounts[j,2] <- nrow(resultstablesubset)
      speciescounts[j,3] <- mean(resultstablesubset$pident)
      speciescounts[j,4] <- max(resultstablesubset$pident)
      speciescounts[j,5] <- min(resultstablesubset$pident)
      speciescounts[j,6] <- resultstablesubset$Sample[1]
      
    }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>3) {
      speciescounts <- speciescounts[1:3,]
      
      
    }
    Euklistdiamondraws[[i]] <- speciescounts
    
  }
  
  if (nrow(ResultstableEuk)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    Euklistdiamondraws[[i]] <- speciescounts
  }
  
  
  
  
  Resultstablebacteria <- subset(Resultstable,Resultstable$superkingdom=="Bacteria")
  
  
  uniquefamilies <- unique(Resultstablebacteria$family)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquefamilies), ncol=6))
  colnames(speciescounts) <- c("Family", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(Resultstablebacteria)>=1) { 
    for ( j in c(1:length(uniquefamilies))) {
      
      resultstablesubset <- subset(Resultstablebacteria,Resultstablebacteria$family==uniquefamilies[j])
      
      speciescounts[j,1] <- resultstablesubset$family[1]
      speciescounts[j,2] <- nrow(resultstablesubset)
      speciescounts[j,3] <- mean(resultstablesubset$pident)
      speciescounts[j,4] <- max(resultstablesubset$pident)
      speciescounts[j,5] <- min(resultstablesubset$pident)
      speciescounts[j,6] <- resultstablesubset$Sample[1]
      
    }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>3) {
      speciescounts <- speciescounts[1:3,]
    }
    bactlistdiamondraws[[i]] <- speciescounts
    
  }
  
  if (nrow(Resultstablebacteria)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    bactlistdiamondraws[[i]] <- speciescounts
  }
  
  ResultstableVirus <- subset(Resultstable,Resultstable$superkingdom=="Viruses")
  
  uniquefamilies <- unique(ResultstableVirus$family)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquefamilies), ncol=6))
  colnames(speciescounts) <- c("Family", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(ResultstableVirus)>=1) { 
    for ( j in c(1:length(uniquefamilies))) {
      
      resultstablesubset <- subset(ResultstableVirus,ResultstableVirus$family==uniquefamilies[j])
      
      speciescounts[j,1] <- resultstablesubset$family[1]
      speciescounts[j,2] <- nrow(resultstablesubset)
      speciescounts[j,3] <- mean(resultstablesubset$pident)
      speciescounts[j,4] <- max(resultstablesubset$pident)
      speciescounts[j,5] <- min(resultstablesubset$pident)
      speciescounts[j,6] <- resultstablesubset$Sample[1]
      
    }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>3) {
      speciescounts <- speciescounts[1:3,]
    }
    virlistdiamondraws[[i]] <- speciescounts
  }
  
  if (nrow(ResultstableVirus)==0) {
    
    speciescounts[1,1]  <- ""
    speciescounts[1,2]  <- 0
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- sampname
    
    virlistdiamondraws[[i]] <- speciescounts
  }
  
  
  
}
# Combine lists
Virusesdiamond <- bind_rows(virlistdiamondraws)
Bacteriadiamond <- bind_rows(bactlistdiamondraws)
Eukaryotesdiamond <- bind_rows(Euklistdiamondraws)


# ggplot graphs 


plot <- ggplot(Virusesdiamond, aes(x=Sample, y=`Reads assigned`, fill=Family)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top three virus families reads/contigs \n  were assigned to, per sample") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7),axis.text.x = element_text(size=3.5)) +
  theme(legend.title = element_text(size=6), legend.text = element_text(size=3.4)) +
  theme(legend.key.size = unit(0.28,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=4.5)) +
  theme(plot.title =  element_text(size=9,hjust = 0.4)) +
  coord_flip() +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot


ggsave(paste0(outtablespath,"Summary_barplot_contigs_virus_families_raw_reads_diamond.pdf"), dpi = 300)

plot

ggsave(paste0(outtablespath,"Summary_barplot_contigs_virus_families_raw_reads_diamond.png"), dpi = 300)


plot <- ggplot(Bacteriadiamond, aes(x=Sample, y=`Reads assigned`, fill=Family)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top 3 Bacterial families reads/contigs \n  were assigned to, per sample") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7),axis.text.x = element_text(size=3.5)) +
  theme(legend.title = element_text(size=6), legend.text = element_text(size=3.4)) +
  theme(legend.key.size = unit(0.28,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=4.5)) +
  theme(plot.title =  element_text(size=9,hjust = 0.4)) +
  coord_flip() +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot


ggsave(paste0(outtablespath,"Summary_barplot_contigs_bacteria_families_raw_reads_diamond.pdf"), dpi = 300)

plot

ggsave(paste0(outtablespath,"Summary_barplot_contigs_bacteria_families_raw_reads_diamond.png"), dpi = 300)


plot <- ggplot(Eukaryotesdiamond, aes(x=Sample, y=`Reads assigned`, fill=Family)) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  #scale_fill_brewer(palette = "Set3") +
  #scale_fill_viridis(discrete = T,name= "Reads assigned as:") +
  ggtitle("Top 3 Eukaryote families reads/contigs \n  were assigned to, per sample") +
  #theme_ipsum() +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7)) +
  theme(axis.title.x = element_text(hjust = 0.5,size = 7),axis.text.x = element_text(size=3.5)) +
  theme(legend.title = element_text(size=6), legend.text = element_text(size=3.4)) +
  theme(legend.key.size = unit(0.28,"line")) +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=4.5)) +
  theme(plot.title =  element_text(size=9,hjust = 0.4)) +
  coord_flip() +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot


ggsave(paste0(outtablespath,"Summary_barplot_contigs_virus_Eukaryote_raw_reads_diamond.pdf"), dpi = 300)


plot

ggsave(paste0(outtablespath,"Summary_barplot_contigs_virus_Eukaryote_raw_reads_diamond.png"), dpi = 300)

}



# Eukaryotes top 20 returned hits sumamry table

sampleEuksummaryfiles <- list.files(path = intablespath, pattern = "top20Eukaryotehits", all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE)


AllsampleEukssummaryfiles  <- list()

for (i in c(1:length(sampleEuksummaryfiles))) {
  
  Resultstable <- read.table(sampleEuksummaryfiles[i],header = TRUE,sep = "\t")
  sample<- Resultstable$Sample
  
  AllsampleEukssummaryfiles[[i]] <- Resultstable
  
}

AllsampleEukssummaryfilesdf <- bind_rows(AllsampleEukssummaryfiles)


AllsampleEukssummaryfilesdf10top <- AllsampleEukssummaryfilesdf %>% 
  filter(species %in% Eukaryotes_top10$Species)

aggtable <- aggregate((AllsampleEukssummaryfilesdf10top$Frequency),by=list(AllsampleEukssummaryfilesdf10top$species,AllsampleEukssummaryfilesdf10top$Sample),sum)
aggtable2table <- aggregate((AllsampleEukssummaryfilesdf10top$average_percent_ident),by=list(AllsampleEukssummaryfilesdf10top$species),mean)
aggtable3table <- aggregate((AllsampleEukssummaryfilesdf10top$length),by=list(AllsampleEukssummaryfilesdf10top$species),mean)


pairwise_matrixEuk <- acast(aggtable, Group.1 ~ Group.2, value.var = "x", fill = 0)

pairwise_matrixEuk <- as.data.frame(pairwise_matrixEuk)

pairwise_matrixEuk$avg_pairwiseidentity <- aggtable2table$x
pairwise_matrixEuk$avg_aligned_length <- aggtable3table$x



write.table(pairwise_matrixEuk,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Eukaryote_species.txt")),sep="\t",row.names=rownames(pairwise_matrixEuk),col.names=colnames(pairwise_matrixEuk))


# Now bacteria 


sampleBacsummaryfiles <- list.files(path = intablespath, pattern = "top20bacterialhits", all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE)


AllsampleBacsummaryfiles <- list()

for (i in c(1:length(sampleBacsummaryfiles))) {
  
  Resultstable <- read.table(sampleBacsummaryfiles[i],header = TRUE,sep = "\t")
  sample<- Resultstable$Sample
  
  AllsampleBacsummaryfiles[[i]] <- Resultstable
  
}

AllsampleBacsummaryfilesdf <- bind_rows(AllsampleBacsummaryfiles)


AllsampleBacsummaryfilesdf10top <- AllsampleBacsummaryfilesdf %>% 
  filter(species %in% Bacteria_top10$Species)



aggtable <- aggregate((AllsampleBacsummaryfilesdf10top$Frequency),by=list(AllsampleBacsummaryfilesdf10top$species,AllsampleBacsummaryfilesdf10top$Sample),sum)
aggtable2table <- aggregate((AllsampleBacsummaryfilesdf10top$average_percent_ident),by=list(AllsampleBacsummaryfilesdf10top$species),mean)
aggtable3table <- aggregate((AllsampleBacsummaryfilesdf10top$length),by=list(AllsampleBacsummaryfilesdf10top$species),mean)



pairwise_matrixbact <- acast(aggtable, Group.1 ~ Group.2, value.var = "x", fill = 0)

pairwise_matrixbact<- as.data.frame(pairwise_matrixbact)

pairwise_matrixbact$avg_pairwiseident <- aggtable2table$x
pairwise_matrixbact$avg_aligned_length <- aggtable3table$x



write.table(pairwise_matrixbact,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Bacterial_species.txt")),sep="\t",row.names=rownames(pairwise_matrixbact),col.names=colnames(pairwise_matrixbact))




# Now Viral 


sampleVirsummaryfiles <- list.files(path = intablespath, pattern = "top20Viralhits_contigs.txt", all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE)


AllsampleVirsummaryfiles <- list()

for (i in c(1:length(sampleVirsummaryfiles))) {
  
  Resultstable <- read.table(sampleVirsummaryfiles[i],header = TRUE,sep = "\t")
  sample<- Resultstable$Sample
  
  AllsampleVirsummaryfiles[[i]] <- Resultstable
  
}

AllsampleVirsummaryfilesdf <- bind_rows(AllsampleVirsummaryfiles)



AllsampleVirsummaryfilesdf10top <- AllsampleVirsummaryfilesdf %>% 
  filter(species %in% Viruses_top10$Species)



aggtable <- aggregate((AllsampleVirsummaryfilesdf10top$Frequency),by=list(AllsampleVirsummaryfilesdf10top$species,AllsampleVirsummaryfilesdf10top$Sample),sum)
aggtable2table <- aggregate((AllsampleVirsummaryfilesdf10top$average_percent_ident),by=list(AllsampleVirsummaryfilesdf10top$species),mean)
aggtable3table <- aggregate((AllsampleVirsummaryfilesdf10top$length),by=list(AllsampleVirsummaryfilesdf10top$species),mean)




pairwise_matrixVir <- acast(aggtable, Group.1 ~ Group.2, value.var = "x", fill = 0)


pairwise_matrixVir <- as.data.frame(pairwise_matrixVir)

pairwise_matrixVir$avg_pairwiseident <- aggtable2table$x
pairwise_matrixVir$avg_aligned_length <- aggtable3table$x



write.table(pairwise_matrixVir,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Viral_species.txt")),sep="\t",row.names=rownames(pairwise_matrixVir),col.names=colnames(pairwise_matrixVir))


write.table(Viruses_top10,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Viral_species_long_table_format.txt")),sep="\t",row.names=FALSE,col.names=colnames(Viruses_top10))

write.table(Eukaryotes_top10,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Eukaryotes_species_long_table_format.txt")),sep="\t",row.names=FALSE,col.names=colnames(Eukaryotes_top10))

write.table(Bacteria_top10,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Bacteria_species_long_table_format.txt")),sep="\t",row.names=FALSE,col.names=colnames(Bacteria_top10))




save.image(paste0(outtablespath,"gather_summary_files_R_environment.Rdata"))
