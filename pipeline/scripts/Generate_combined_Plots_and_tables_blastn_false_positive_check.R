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

# Now bring in the Contig assigned read files to generate figures for top returned species. 

# Binary yes no to whether to read in the no host freq summary or the full freq summary

# It can read on empty files, but then there is no sample name which was stored in the files otherwise

# Top 5 of each superkingdom as separate lists per sample. run files to create and assign each. Then rbind list values into combined data frame for 
# combined figure generation.


  
  summaryreturnedhitsfromcontigs <- list.files(path = intablespath, pattern = "summarycontighits_assigned_assembly_including_blastn_false_positive_check", all.files = FALSE,
                                               full.names = TRUE, recursive = FALSE,
                                               ignore.case = FALSE, include.dirs = FALSE)
  


Euklist <- list()
bactlist <- list()
virlist <- list()



for (i in c(1:length(summaryreturnedhitsfromcontigs))) {
  
  Resultstable <- read.table(summaryreturnedhitsfromcontigs[i],header = TRUE,sep = "\t")
  sampname <- gsub(".*\\/" , x=summaryreturnedhitsfromcontigs[i],replacement="")
  sampname <- gsub("_summarycontig.*",x=sampname,replacement="")
  Resultstable <- subset(Resultstable, !(is.na(Resultstable$species)))
  # Need to remove all NA rows from data
  
  ResultstableEuk <- subset(Resultstable,Resultstable$superkingdom=="Eukaryota")
  
  uniquespecies <- unique(ResultstableEuk$subspecies)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=7))
  colnames(speciescounts) <- c("Species","subspecies", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(ResultstableEuk)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(ResultstableEuk,ResultstableEuk$subspecies==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- resultstablesubset$subspecies[1]
    speciescounts[j,3] <- sum(resultstablesubset$freq)
    speciescounts[j,4] <- mean(resultstablesubset$percentident)
    speciescounts[j,5] <- max(resultstablesubset$percentident)
    speciescounts[j,6] <- min(resultstablesubset$percentident)
    speciescounts[j,7] <- sampname    
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }
  Euklist[[i]] <- speciescounts
  
}
  
  if (nrow(ResultstableEuk)==0) {
    
    speciescounts[1,1]  <- "NA"
    speciescounts[1,2]  <- "NA"
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- 0
    speciescounts[1,7]  <- sampname
    
    Euklist[[i]] <- speciescounts
  }
  
  
  
  
  Resultstablebacteria <- subset(Resultstable,Resultstable$superkingdom=="Bacteria")
  
  
  uniquespecies <- unique(Resultstablebacteria$subspecies)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=7))
  colnames(speciescounts) <- c("Species", "subspecies", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Sample")
  
  if (nrow(Resultstablebacteria)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(Resultstablebacteria,Resultstablebacteria$subspecies==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- resultstablesubset$subspecies[1]
    speciescounts[j,3] <- sum(resultstablesubset$freq)
    speciescounts[j,4] <- mean(resultstablesubset$percentident)
    speciescounts[j,5] <- max(resultstablesubset$percentident)
    speciescounts[j,6] <- min(resultstablesubset$percentident)
    speciescounts[j,7] <- sampname
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }
  bactlist[[i]] <- speciescounts
  
}
  
  if (nrow(Resultstablebacteria)==0) {
    
    speciescounts[1,1]  <- "NA"
    speciescounts[1,2]  <- "NA"
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- 0
    speciescounts[1,7]  <- sampname
    
    bactlist[[i]] <- speciescounts
  }
  


  ResultstableVirus <- subset(Resultstable,Resultstable$superkingdom=="Viruses")
  
  uniquespecies <- unique(ResultstableVirus$subspecies)
  speciescounts <- as.data.frame(matrix(nrow = length(uniquespecies), ncol=12))
  colnames(speciescounts) <- c("Species", "subspecies", "Reads assigned", "mean percent identity to hit","max percent identity to hit","min percent identity to hit","Blastn alternate kingdom identified","Blastn alternate species identified","Blastn alternate subspecies identified","Blastn alternate hit average pairwise identity","Blastn alternate hit average alignment length","Sample")
  
  if (nrow(ResultstableVirus)>=1) { 
  for ( j in c(1:length(uniquespecies))) {
    
    resultstablesubset <- subset(ResultstableVirus,ResultstableVirus$subspecies==uniquespecies[j])
    
    speciescounts[j,1] <- resultstablesubset$species[1]
    speciescounts[j,2] <- resultstablesubset$subspecies[1]
    speciescounts[j,3] <- sum(resultstablesubset$freq)
    speciescounts[j,4] <- mean(resultstablesubset$percentident)
    speciescounts[j,5] <- max(resultstablesubset$percentident)
    speciescounts[j,6] <- min(resultstablesubset$percentident)


    if (sum(resultstablesubset$blastn_alternate_superkingdom_id != "Viruses", na.rm = TRUE) >=1) {
	resultstablesubsetalts <- subset(resultstablesubset,resultstablesubset$blastn_alternate_superkingdom_id != "Viruses")

    speciescounts[j,7] <- resultstablesubsetalts$blastn_alternate_superkingdom[1]
    speciescounts[j,8] <- resultstablesubsetalts$blastn_alternate_species[1]
    speciescounts[j,9] <- resultstablesubsetalts$blastn_alternate_subspecies[1]


	}
    if (sum(resultstablesubset$blastn_alternate_superkingdom != "Viruses", na.rm = TRUE) <1) {

		if (sum(resultstablesubset$blastn_false_positive_check =="Yes",na.rm=TRUE)>=1) {

		    speciescounts[j,7] <- "Viral species returned from Blastn"
		    speciescounts[j,8] <- resultstablesubset$blastn_alternate_species[1]
		    speciescounts[j,9] <- resultstablesubset$blastn_alternate_subspecies[1]
		}

		if (sum(resultstablesubset$blastn_false_positive_check =="Yes",na.rm=TRUE)<1) {

		    speciescounts[j,7] <- "Blastn returned no species"
		    speciescounts[j,8] <- "None"
		    speciescounts[j,9] <- "None"
		}

	}

    speciescounts[j,10] <- mean(as.numeric(resultstablesubset$`blastn_alternate_percentident`))
    speciescounts[j,11] <- mean(as.numeric(resultstablesubset$`blastn_alternate_alignment_length`))
    speciescounts[j,12] <- sampname
    
  }
    speciescounts <- speciescounts[order(-speciescounts$`Reads assigned`),]
    if(nrow(speciescounts)>10) {
      speciescounts <- speciescounts[1:10,]
    }

  virlist[[i]] <- speciescounts

  }
  
  if (nrow(ResultstableVirus)==0) {
    
    speciescounts[1,1]  <- "NA"
    speciescounts[1,2]  <- "NA"
    speciescounts[1,3]  <- 0
    speciescounts[1,4]  <- 0
    speciescounts[1,5]  <- 0
    speciescounts[1,6]  <- 0
    speciescounts[1,7]  <- "None"
    speciescounts[1,8]  <- "None"
    speciescounts[1,9]  <- "None"
    speciescounts[1,10]  <- 0
    speciescounts[1,11]  <- 0
    speciescounts[1,12] <- sampname

    
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
Viruses_sums <- aggregate((Viruses$`Reads assigned`),by=list(Viruses$subspecies),sum)
Viruses_sums <- Viruses_sums[order(-Viruses_sums$`x`),]

if (nrow(Viruses_sums)>=10) {
  Viruses_sumstop10 <- Viruses_sums[1:10,]
}
if (nrow(Viruses_sums)<10) {
  Viruses_sumstop10 <- Viruses_sums
}

Viruses_top10 <- Viruses %>% 
  filter(subspecies %in% Viruses_sumstop10$Group.1)

Viruses_sumstop10$Group.1


# Bacteria
Bacteria_sums <- Bacteria[order(-Bacteria$`Reads assigned`),]

Bacteria_sums <- aggregate((Bacteria$`Reads assigned`),by=list(Bacteria$subspecies),sum)
Bacteria_sums <- Bacteria_sums[order(-Bacteria_sums$`x`),]

if (nrow(Bacteria_sums)>=10) {
  Bacteria_sumstop10 <- Bacteria_sums[1:10,]
}
if (nrow(Bacteria_sums)<10) {
  Bacteria_sumstop10 <- Bacteria_sums
}

Bacteria_top10 <- Bacteria %>% 
  filter(Species %in% Bacteria_sumstop10$Group.1)


# Eukaryotes
Eukaryotes_sums <- Eukaryotes[order(-Eukaryotes$`Reads assigned`),]

Eukaryotes_sums <- aggregate((Eukaryotes$`Reads assigned`),by=list(Eukaryotes$subspecies),sum)
Eukaryotes_sums <- Eukaryotes_sums[order(-Eukaryotes_sums$`x`),]

if (nrow(Eukaryotes_sums)>=10) {
  Eukaryotes_sumstop10 <- Eukaryotes_sums[1:10,]
}
if (nrow(Eukaryotes_sums)<10) {
  Eukaryotes_sumstop10 <- Eukaryotes_sums
}

Eukaryotes_top10 <- Eukaryotes %>% 
  filter(Species %in% Eukaryotes_sumstop10$Group.1)






# ggplot graphs 

Viruses_top10$`Reads assigned` <- as.numeric(as.character(Viruses_top10$`Reads assigned`))
Bacteria_top10$`Reads assigned` <- as.numeric(as.character(Bacteria_top10$`Reads assigned`))
Eukaryotes_top10$`Reads assigned` <- as.numeric(as.character(Eukaryotes_top10$`Reads assigned`))


Viruses_top10$Species2 <- NA

for (a in c(1:nrow(Viruses_top10))) {


if (Viruses_top10$`Blastn alternate kingdom identified`[a] !="Viral species returned from Blastn") {

Viruses_top10$Species2[a] <- paste0(Viruses_top10$subspecies[a],"*")

}

else {
Viruses_top10$Species2[a] <- Viruses_top10$subspecies[a]
}

}


# ggplot graphs 


plot <- ggplot(Viruses_top10, aes(x=Sample, y=`Reads assigned`, fill=Species2)) +
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
  labs(fill = "Species") +
  theme(axis.title.y = element_blank(), axis.text.y =  element_text(size=16.5)) +
  theme(plot.title =  element_text(size=19,hjust = 0.4)) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(labels = scales::number_format()) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

plot

ggsave(paste0(outtablespath,"Summary_barplot_top_virus_species_from_contigs_including_false_positive_check.pdf"), dpi = 330)

plot 

ggsave(paste0(outtablespath,"Summary_barplot_top_virus_species_from_contigs_including_false_positive_check.png"), dpi = 330)


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
  filter(subspecies %in% Viruses_top10$subspecies)




aggtable <- aggregate((AllsampleVirsummaryfilesdf10top$Frequency),by=list(AllsampleVirsummaryfilesdf10top$subspecies,AllsampleVirsummaryfilesdf10top$Sample),sum)
aggtable2table <- aggregate((AllsampleVirsummaryfilesdf10top$average_percent_ident),by=list(AllsampleVirsummaryfilesdf10top$subspecies),mean)
aggtable3table <- aggregate((AllsampleVirsummaryfilesdf10top$length),by=list(AllsampleVirsummaryfilesdf10top$subspecies),mean)



pairwise_matrixVir <- acast(aggtable, Group.1 ~ Group.2, value.var = "x", fill = 0)


pairwise_matrixVir <- as.data.frame(pairwise_matrixVir)

pairwise_matrixVir$avg_pairwiseident <- aggtable2table$x
pairwise_matrixVir$avg_aligned_length <- aggtable3table$x







Virus_final_mixed_table <- as.data.frame(matrix(nrow=nrow(pairwise_matrixVir), ncol=(ncol(pairwise_matrixVir))))

colnames(Virus_final_mixed_table) <- colnames(pairwise_matrixVir) 
rownames(Virus_final_mixed_table) <- rownames(pairwise_matrixVir)

for ( i in c(1:(ncol(Virus_final_mixed_table)-2))) {

grep(colnames(Virus_final_mixed_table[i]),Viruses_top10$Sample) -> index1 

Viruses_top10sampleX <-Viruses_top10[index1,]

for (j in c(1:nrow(Virus_final_mixed_table))) {


rownames(Virus_final_mixed_table)[j] -> Speciesnametest

Viruses_top10sampleXvirsubset <- subset(Viruses_top10sampleX,Viruses_top10sampleX$subspecies==Speciesnametest)

if (nrow(Viruses_top10sampleXvirsubset) >=1) {

if (Speciesnametest != Viruses_top10sampleXvirsubset$Species2) {

Virus_final_mixed_table[j,i] <- paste0(pairwise_matrixVir[j,i],"*")

}

if (Speciesnametest == Viruses_top10sampleXvirsubset$Species2) {

Virus_final_mixed_table[j,i] <- pairwise_matrixVir[j,i]

}
}

}

j=1

}

Virus_final_mixed_table$avg_pairwiseident <- aggtable2table$x
Virus_final_mixed_table$avg_aligned_length <- aggtable3table$x

Virus_final_mixed_table[is.na(Virus_final_mixed_table)] <- 0



write.table(Virus_final_mixed_table,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Viral_species_blastn_false_positive_check.txt")),sep="\t",row.names=rownames(pairwise_matrixVir),col.names=colnames(pairwise_matrixVir))

Viruses_top10red <- Viruses_top10[1:10,]





write.table(Viruses_top10red,file=(paste0(outtablespath,"Combined_samples_top_hits_to_Viral_species_long_table_format_blastn_false_positive_check.txt")),sep="\t",row.names=FALSE,col.names=colnames(Viruses_top10red))



save.image(paste0(outtablespath,"gather_summary_files_R_environment_Blastn_false_positive_check.Rdata"))
