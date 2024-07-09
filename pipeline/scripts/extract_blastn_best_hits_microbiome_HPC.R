# Script for subsetting and extracting valuble information from Blastn microbiomes
# Not optimised for the microbiome. Microbiome only returns top hit so the multi species sections are unnecessary. Largely copy from original blastn analysis script 

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
parser$add_argument('--output', '-o', help= 'Output file for assigned contigs')
parser$add_argument('--abundances', '-a', help= 'Output file path for abundance hit tables')
parser$add_argument('--programdir ', '-p', help= 'working program directory')
parser$add_argument('--savdir ', '-s', help= 'working directory for saving contig matches')
parser$add_argument('--savcontig ', '-S', help= 'filename for saving contig matches')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--threads', '-t', help= 'Number of threads')
parser$add_argument('--Acc', '-A', help= 'Accessiontaxa filepath')
parser$add_argument('--Log', '-l', help= 'Log of data')
parser$add_argument('--Accnode', '-N', help= 'Accessiontaxa name_node filepath')


xargs<- parser$parse_args()
blastn_output <- read.table(xargs$inputblastn, header = FALSE, sep = "\t",fill=TRUE,comment.char="")
outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$abundances
outtablespath <- paste0(outtablespath1,outtablespath2)
AccessionNamenode <- xargs$Accnode

NAMES <- xargs$name
n.cores <- xargs$threads
Accessionfile <- xargs$Acc

assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
assigned_contigs <- paste0(assigned_contigs1,assigned_contigs2)

contigsbasedir <- xargs$savdir
contigsvdir <- paste0(assigned_contigs1,contigsbasedir)

#contigsvdir <-xargs$savdir
svcontigname <-xargs$savcontig

log <- xargs$Log


 

# Set to correct directory.
# sink all cats and prints to the snakemake log files 
 sink(log)

if (!file.exists(outtablespath)){

dir.create(outtablespath, recursive = FALSE, mode = "0777")
}

cat(paste0("printing results tables to", outtablespath))
cat(paste0("printing matched contigs to", assigned_contigs))

  # give proper column names
  colnames(blastn_output) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "staxids", "stitle", "qcovhsp")
  
  cat(paste0("Succesfully read in Diamond file for sample", NAMES,"\n"))
  
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
  contigsassigned <- as.data.frame(matrix(nrow=length(contigs),ncol = ncol(blastn_output)))
  colnames(contigsassigned) <- colnames(blastn_output)
  
  # Loop serves the functions of 
  # 1. subsettting data to each contig (all X hits returned are subset into new df for each contig) 
  # 2. determining whether  the hits returned similar in bitscore or different
  # 3. if similar, determine wha is likely the best species based on the frequency of hits to each
  # taxid
  
a <- Sys.time()  

  for (i in (c(1:length(contigs)))) {
    
    workingcontigA <- subset(blastn_output,blastn_output$qseqid==contigs[i]) 
    workingcontig <- workingcontigA[order(-workingcontigA$bitscore),]
    namecounter=1
    # If here required in case Diamond blastx only returns a single hit 
    if (nrow(workingcontig)==1 ) {
      
      contigsassigned[i,] <- workingcontig[1,]
    }
    
    # else if here tests to see if the top hit is more than 2% higher bitscore than the next hit
    # If it is, only the top hit information is taken. The 2% is arbitrary and should look further 
    # into other thresholds. The idea is to compare whether multiple species are closely aligned vs
    # One species hit returned is much closer
    else if (workingcontig$bitscore[1] > 1.02*(workingcontig$bitscore[2]) )  {
      
      contigsassigned[i,] <- workingcontig[1,]
      
    }
    # else if here is for when the bitscores are  <=2% different. This compares
    # all of the different hits to the top hit and collects the taxids of all that are <5% different
    # It tallys the total taxids for the contig (e.g., 5 hits are mosquito X, 2 are mosquito Y).
    # It then extracts the highest bitscore hit for the most frequently occuring taxid e.g. Mosquito X
    
    else if (workingcontig$bitscore[1] <= 1.02*(workingcontig$bitscore[2])) {
      
      
      
      for (j in (c(1:nrow(workingcontig)))) {
        
        
        
        if (workingcontig$bitscore[1] <=1.02*workingcontig$bitscore[j]) {
          
          
          namesorter[namecounter] <- workingcontig$staxids[j]
          
          
          namecounter=namecounter+1
          
          
          
        }
      }
      namesorter <- namesorter[1:(namecounter-1)]
      namehits <- plyr::count(namesorter)
      
      namehitsidx <- order(namehits$freq,decreasing = TRUE,na.last = TRUE)
      
      
      namehittaxid <- namehits[namehitsidx[1],1]
      
      stopval=0
      
      if(length(namehitsidx)==1) {
        
        contigsassigned[i,] <- workingcontig[1,]
        
        
      }
      
      else if (length(namehitsidx)>=2) {
        
        for (b in (c(1:nrow(workingcontig)))) {
          
          if (namehits$freq[namehitsidx[1]]==namehits$freq[namehitsidx[2]]) {
            
            contigsassigned[i,] <- workingcontig[1,]
            
            
          }
          
          else if (workingcontig$staxids[b]==namehittaxid & stopval==0 & namehits$freq[namehitsidx[1]]!=namehits$freq[namehitsidx[2]]) {
            
            contigsassigned[i,] <- workingcontig[b,]
            
            
            stopval=1
            
          }
          
          
          
        }
      }
    }
    
    
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

taxids <- as.data.frame(matrix(nrow=nrow(contigsassigned),ncol=8))
taxidsunique <- as.data.frame(matrix(nrow=nrow(contigsassignedunique),ncol=8))
taxids[,1] <-contigsassigned$staxidreduced
taxidsunique[,1] <-contigsassignedunique$staxidreduced

taxidsunique[,2:8] <- foreach(i = 1:nrow(contigsassignedunique), .combine = 'rbind') %dopar% {

taxidsunique[i,2:8] <- taxonomizr::getTaxonomy(contigsassignedunique$staxidreduced[i], sqlFile=AccessionNamenode)

}


d <- Sys.time()

	cat(paste0(NAMES," Finished taxid conversion to taxonomy ", "\n"))
	cat(paste0(Sys.time(), "\n"))
	cat(paste0(" Time taken: ", difftime(d,c), "\n"))




taxidsunique<- as.data.frame(taxidsunique, ncol=8)


e <- Sys.time()

taxids$V1 <- gsub(pattern = ".",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ".",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = "(",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = "(",replacement = "",x = taxidsunique$V1,fixed = TRUE)
taxids$V1 <- gsub(pattern = ")",replacement = "",x = taxids$V1,fixed = TRUE)
taxidsunique$V1 <- gsub(pattern = ")",replacement = "",x = taxidsunique$V1,fixed = TRUE)

colnames(taxidsunique) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species") 

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
colnames(taxids) =c("staxidreduced", "superkingdom", "phylum", "class", "order", "family", "genus", "species") 

contigsassigned$superkingdom <- taxids$superkingdom
contigsassigned$phylum<- taxids$phylum
contigsassigned$class<- taxids$class
contigsassigned$order<- taxids$order
contigsassigned$family<- taxids$family
contigsassigned$genus<- taxids$genus
contigsassigned$species<- taxids$species

contigsassigned$stitle = substr(contigsassigned$stitle,1,50)


      cat(paste0(NAMES," Finished reindexing of taxa results ", "\n"))
      cat(paste0(Sys.time(), "\n"))

  
# generate summary stats for blast
  
  # Now this section is used to update the original assembled contig names so that they include
  # the additional species information (This step may be superfluous later on but it may also be useful
  # to compare to see how genome binning works compared to direct assignment)
  
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



superkingdom <- plyr::count(contigsassigned$superkingdom)



  write.table(Eukaryotessinglespspcounts, file=(paste0(outtablespath,"/",NAMES,"_Eukaryotes_sp_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Eukaryotessinglespspcounts),sep="\t")
  write.table(Bacteriasinglespspcounts, file=(paste0(outtablespath,"/",NAMES,"_Bacteria_species_hits.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(Bacteriasinglespspcounts),sep="\t")
  write.table(contigsassigned, file=(paste0(outtablespath,"/",NAMES,"_Contigsallinformationassignment.txt")), quote=FALSE, row.names=FALSE, col.names=colnames(contigsassigned),sep="\t")


cat(paste0("Finished individual ",NAMES,"\n"))

