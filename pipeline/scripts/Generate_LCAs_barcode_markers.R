library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(htmltools)
library(htmlwidgets)
library(magrittr)
library(argparse)
library(hrbrthemes)
library(sankeyD3)
library(pavian)
library("devtools")


if (!require("d3treeR", character.only = TRUE)) {
  devtools::install_github("timelyportfolio/d3treeR",upgrade="never",dependencies=FALSE)
  library("d3treeR")
} else {
  library("d3treeR")
}




parser <- ArgumentParser(description= 'Summarising results filtering and assembly')

parser$add_argument('--programdir', '-a', help= 'Input taxonomy file')
parser$add_argument('--outputpath', '-b', help= 'Input taxonomy file')
parser$add_argument('--inputtax', '-g', help= 'Input taxonomy file')
parser$add_argument('--name', '-n', help= 'Name of sample')
parser$add_argument('--host', '-c', help= 'Whether one or multipe host species are expected')
parser$add_argument('--marker', '-p', help= 'Marker being analysed')

xargs<- parser$parse_args()

inputtaxa <- xargs$inputtax

NAMES <- xargs$name

basepath <- xargs$programdir
resultspath <- xargs$outputpath

outtablespath <- paste0(basepath,resultspath)

marker <- xargs$marker

Hostdetect <- xargs$host

myfile <- pavian::read_report(xargs$inputtax)


# sankey function

build_sankey_network <- function(my_report, taxRanks =  c("D","K","P","C","O","F","G","S"), maxn=10,
                                 zoom = F, title = NULL,
                                 ...) {
  stopifnot("taxRank" %in% colnames(my_report))
  if (!any(taxRanks %in% my_report$taxRank)) {
    warning("report does not contain any of the taxRanks - skipping it")
    return()
  }
  my_report <- subset(my_report, taxRank %in% taxRanks)
  my_report <- plyr::ddply(my_report, "taxRank", function(x) x[utils::tail(order(x$cladeReads,-x$depth), n=maxn), , drop = FALSE])
  
  my_report <- my_report[, c("name","taxLineage","taxonReads", "cladeReads","depth", "taxRank")]
  
  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  #my_report$name <- sub("^-_root.", "", my_report$name)
  
  splits <- strsplit(my_report$taxLineage, "\\|")
  
  ## for the root nodes, we'll have to add an 'other' link to account for all cladeReads
  root_nodes <- sapply(splits[sapply(splits, length) ==2], function(x) x[2])
  
  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]
  
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n=2))), stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  links$value <- my_report[sel,"cladeReads"]
  
  my_taxRanks <- taxRanks[taxRanks %in% my_report$taxRank]
  taxRank_to_depth <- stats::setNames(seq_along(my_taxRanks)-1, my_taxRanks)
  
  
  nodes <- data.frame(name=my_report$name,
                      depth=taxRank_to_depth[my_report$taxRank],
                      value=my_report$cladeReads,
                      stringsAsFactors=FALSE)
  
  for (node_name in root_nodes) {
    diff_sum_vs_all <- my_report[my_report$name == node_name, "cladeReads"] - sum(links$value[links$source == node_name])
    if (diff_sum_vs_all > 0) {
      nname <- paste("other", sub("^._","",node_name))
      #nname <- node_name
      #links <- rbind(links, data.frame(source=node_name, target=nname, value=diff_sum_vs_all, stringsAsFactors = FALSE))
      #nodes <- rbind(nodes, nname)
    }
  }
  
  names_id = stats::setNames(seq_len(nrow(nodes)) - 1, nodes[,1])
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[links$source != links$target, ]
  
  nodes$name <- sub("^._","", nodes$name)
  links$source_name <- nodes$name[links$source + 1]
  
  if (!is.null(links))
    sankeyD3::sankeyNetwork(
      Links = links,
      Nodes = nodes,
      doubleclickTogglesChildren = TRUE,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      NodeGroup = "name",
      NodePosX = "depth",
      NodeValue = "value",
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 12,
      nodeCornerRadius = 5,
      units = "cladeReads",
      fontSize = 12,
      iterations = maxn * 100,
      align = "none",
      highlightChildLinks = TRUE,
      orderByPath = TRUE,
      scaleNodeBreadthsByString = TRUE,
      zoom = zoom,
      ...
    )
}



#

#


#


#

#

# Start of analysis code

myfile <- pavian::read_report(inputtaxa)


if (!(is.null(myfile))) {
  
  if ( (sum(grepl("D",x = myfile$taxRank)) >=2) | sum(grepl("K",x = myfile$taxRank)) >=1) {
    
    print(c(1:10) )
    
    
    sankeynetwork <- build_sankey_network(myfile,taxRanks = c("D","K","F","G","S"),zoom = FALSE)
    
    paste0("saving output figure to the following ", outtablespath,NAMES,"_taxonomic_classifications", marker, ".html")
    
    
    saveNetwork <- function(network, file, selfcontained = TRUE) {
      htmlwidgets::saveWidget(network, file, selfcontained)
    }
    
    
    
    sankeynetwork %>% saveNetwork(file =paste0(outtablespath,NAMES,"_taxonomic_classifications", marker, ".html"))
  }
}


# Now to identify the most probable host species/blood meal targets

# subset to animals
Animallist <-grep("k_Metazoa",x = myfile$taxLineage)
Animals <- myfile[Animallist,]

# Detect top returned families and order from most reads assigned to fewest
family <- subset(Animals,Animals$taxRank=="F")
speciesrank <- subset(Animals,Animals$taxRank=="S")

if(nrow(speciesrank>=1)) {
  
  if(nrow(family)>=1) {
    
    familyordered <- family[order(-family$cladeReads),]
    
    # Detect top returned family name 
    hostfamily <- familyordered$name[1]
    
    # return all taxonomic hits within the top family
    Famlist <-grep(hostfamily,x = Animals$taxLineage)
    topfamreturns <- Animals[Famlist,]
    
    # Detect top reterned genera within the top returned family and order from most reads assigned to fewest
    genus <- subset(topfamreturns,topfamreturns$taxRank=="G")
    genusordered <- genus[order(-genus$cladeReads),]
    # detect top returned genus name
    hostgenus <- genusordered$name[1]
    
    # Detect top returned species within the top returned family (across all genera)
    species <- subset(topfamreturns,topfamreturns$taxRank=="S")
    speciesordered <- species[order(-species$cladeReads),]
    
    # define the top two genera within the top family separately for comparisons
    Topgenera1 <- genusordered[1,]
    Topgenera2 <- genusordered[2,]
    
    #return all taxonomic hits for the top genus within the top family
    speciesgenera1list <- grep(Topgenera1$name[1],x = Animals$taxLineage)
    topgenerareturns1 <- Animals[speciesgenera1list,]
    # Detect top returned species within the top returned genus
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    #return all taxonomic hits for the second top genus within the top family
    speciesgenera2list <- grep(Topgenera2$name[1],x = Animals$taxLineage)
    topgenerareturns2 <- Animals[speciesgenera2list,]
    
    # Detect top returned species within the second top returned genus
    topspeciesgenera2 <- subset(topgenerareturns2,topgenerareturns2$taxRank=="S")
    topspeciesgenera2ordered <- topspeciesgenera2[order(-topspeciesgenera2$cladeReads),]
    
    # Now need to do pairwise comparisons to make sure one species is being detected sufficiently more frequently to justify their classification as the host species
    
    # Two main subsets here. 1, is there more than one detected fam/genus/sp per subsetted class above. Willr equire a bit of heirarchical analysis here
    if (nrow(familyordered) >=2) {
      
      if (familyordered$cladeReads[1] >= 3*familyordered$cladeReads[2]) {
        
        paste0( " Over 3 times as many reads were assigned to the top family when compared to the next family. This family likely contains the host species ")
        downscalefam <- "no"
        
        paste0( " Now inferring whether there is sufficient resolution to determine the likely genus of the species ")
      }
      
      
      if (familyordered$cladeReads[1] < 3*familyordered$cladeReads[2]) {
        
        paste0( " There is insufficient differentiation between the top returned family and second top returned family. This analysis will return the 
            top hitting species across the whole family ")
        
        downscalefam <- "yes"
      }
    } # Determine whether families are sufficiently differentiated
    
    if (nrow(familyordered) ==1) { 
      
      paste0( " Only one family was detected in this analysis and so is likely the host species family ")
      downscalefam <- "no"
      
    } # quick prep loop for below to just allow for downscale fam to equal no and the below loops to apply for when there is only 1 fam to.
    
    
if (downscalefam == "no") {
  
  if (nrow(genusordered)>=2) {
    
    if (genusordered$cladeReads[1] >= 3*genusordered$cladeReads[2]) {
      
      paste0( " Over 3 times as many reads were assigned to the top genus within the top family when compared to the next genus. 
          This genus likely contains the host species, however the top species across the top genera (assuming species were observed in each genus) will also be compared 
            to find potential discrepencies in case the top hit species is not within the top genus ")
      
      
      
      downscalegenus <- "no"
      downscalegenus2 <- "no"
      
    }
    if (genusordered$cladeReads[1] < 3*genusordered$cladeReads[2]) {
      paste0( " Less than 3 times more reads were assigned to the top genus within the top family when compared to the next genus.
             It may either not be feasible to reliably infer the genus of the host using this marker or it may be that multiple genera were sampled.
             Returning the top species of each id at a minimum of 5% assignment")
      
      downscalegenus <- "yes"
      downscalegenus2 <- "no"
      
    }
    
    
  }
  
  
  
  if (nrow(genusordered)==1) {
    
    paste0( " Only one genus within the top family was detected. Will now determine the top species for this genus and save corresponding stats for
            barcode combination analysis " )
    
    
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    
    
    topspeciesgenera1top5percent <- subset(topspeciesgenera1ordered,topspeciesgenera1ordered$percentage >=1)
    
    if (nrow(topspeciesgenera1top5percent) >=1) {
      
      hostspecies <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=1))
      hostspeciespercents <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=3))
      colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
      
      
      for (i in c(1:nrow(hostspecies))) {
        hostspecies$V1[i] <- gsub(pattern="s_",x=topspeciesgenera1top5percent$name[i],replacement="")
      }
      
      
      hostspeciespercents$species <- hostspecies$V1
      hostspeciespercents$percentage_reads_total_from_marker <- topspeciesgenera1top5percent$percentage
      
      hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
      
      
      write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
      write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
      
      downscalegenus <- "NA"
      downscalegenus2 <- "NA"
      
    }
    
    if (nrow(topspeciesgenera1top5percent) == 0) {
      
      hostspecies <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1ordered),ncol=1))
      hostspeciespercents <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1ordered),ncol=3))
      colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
      
      paste(" No Eukaryote species was found to represent more than 5% of reads. if over 25% of reads are assigned to one genus, LCA marker will be reported \n if not, host identification will be treated as null")
      
      
      if (sum(genusordered$percentage) < 25) {
        
        hostspeciespercents[1,1] <- NA
        hostspeciespercents[1,2] <- 0
        hostspeciespercents[1,3] <- 0
        
      }
      
      for (i in c(1:nrow(hostspecies))) {
        hostspecies$V1[i] <- gsub(pattern="s_",x=topspeciesgenera1ordered$name[i],replacement="")
      }
      
      if (sum(genusordered$percentage) > 25) {
        
        
        hostspeciespercents$species <- hostspecies$V1
        hostspeciespercents$percentage_reads_total_from_marker <- topspeciesgenera1ordered$percentage
        
        hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
        
      }      
      write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
      write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
      
      downscalegenus <- "NA"
      downscalegenus2 <- "NA"
      
    }
    
    
    
    
    
  } # returns top soecies for only top genus. No warnings
  
  
  if (nrow(topspeciesgenera2ordered)>=1) {
    
    if (topspeciesgenera1ordered$percentage[1] < 3* topspeciesgenera2ordered$percentage[1])   {
      
      paste0( " The top species from the second most common genus has more than 33% of the reads assigned as the top species from the top genus. 
            This suggests either the marker utilised fails to reliably resolve species delineation into true monophyletic groups or that multiple species from 
            different genera were sequenced. Returning all top species of each id from the top two genera which meet a minimum of 5% assignment of total reads" )
      
      downscalegenus <- "yes"
      downscalegenus2 <- "yes"
      
      
      
    } # if top species genus 1  doesn't have more hits than top species 
    #genus 2 then prep downstream file writing for host issues
    
    if (topspeciesgenera1ordered$percentage[1] >= 3* topspeciesgenera2ordered$percentage[1])   {
      
      paste0( " Over 3 times as many reads were assigned to the top species within the top family when compared to the top species from the next highest genus. 
          It is likely that the data reliably inferred the correct genus and the marker is sufficiently resolving to inform host genus. Only the top
              species of the top genus will be returned")
      
      downscalegenus <- "no"
      downscalegenus2 <- "no"
      
    } # Confirm top species of top genus has more hits than top species 
    # of second top genus so only one genus needs to be returned
    
    
  } # Analysis step determining downscalegenus and downscalegenus2 variables
  
  
  
  
  
  
  
  
  
  
  
  
  if (nrow(topspeciesgenera2ordered)==0 && downscalegenus != "NA") {
    
    paste0( " The second top returned genus within the top family had not species level classifications. " )
    
    
    
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    
    
    topspeciesgenera1top5percent <- subset(topspeciesgenera1ordered,topspeciesgenera1ordered$percentage >=1)
    
    hostspecies <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=1))
    hostspeciespercents <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=3))
    colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
    
    
    for (i in c(1:nrow(hostspecies))) {
      hostspecies$V1[i] <- gsub(pattern="s_",x=topspeciesgenera1top5percent$name[i],replacement="")
    }
    
    
    hostspeciespercents$species <- hostspecies$V1
    hostspeciespercents$percentage_reads_total_from_marker <- topspeciesgenera1top5percent$percentage
    
    hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
    
    if (downscalegenus=="yes") {
      hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
      hostspecies[1] <- "Can not reliably determine host species for this marker due to divergence at the genus level.
    A top species was observed, and it belonged to the top genus and family, however a large number of reads are assigned to a secondary genus with no
    species classifications suggesting limitations of the resultion of this marker for this genus.
    See the top_host_species_additional_stats for a list of probable species from the top returned family "
      
      
    }
    
    
    write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
    write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
    
    
    downscalegenus <- "NA"
    downscalegenus2 <- "NA"
    
    
  } # if only 1 genus was present and only one family likely then save top species for top 
  
  # Genus to correct files. Possible second pathway for if a secondary genus had high hits but not to a species within that genus (gives warning)
  
  
  if (downscalegenus=="yes" && downscalegenus2=="no") {
    
    
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    
    downscalefam_species <-  subset(topfamreturns,topfamreturns$taxRank=="S")
    
    downscalefam_speciestop5percent <- subset(downscalefam_species,downscalefam_species$percentage >=5)
    
    hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
    hostspeciespercents <- as.data.frame(matrix(nrow=nrow(downscalefam_speciestop5percent),ncol=3))
    colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
    
    hostspecies[1] <- "A likely top family was observed but the two top genera each had similar assignment rates. Can not reliably determine species. See the top_host_species_additional_stats for a list of probable species from the top returned family " 
    
    for (i in c(1:nrow(hostspeciespercents))) {
      hostspeciespercents$species[i] <- gsub(pattern="s_",x=downscalefam_speciestop5percent$name[i],replacement="")
    }
    
    hostspeciespercents$percentage_reads_total_from_marker <- downscalefam_speciestop5percent$percentage
    
    hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
    
    
    write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
    write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
    
    
    downscalegenus <- "NA"
    downscalegenus2 <- "NA"
    
    
  } # returns top species from top family and warning indicating multiple genera were potentially the host genera
  
  
  if (downscalegenus2=="yes") {
    
    
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    
    downscalefam_species <-  subset(topfamreturns,topfamreturns$taxRank=="S")
    
    downscalefam_speciestop5percent <- subset(downscalefam_species,downscalefam_species$percentage >=1)
    
    hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
    hostspeciespercents <- as.data.frame(matrix(nrow=nrow(downscalefam_species),ncol=3))
    colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
    
    hostspecies[1] <- "A likely top genus was observed however a species from the second top genera within the family had similar assignment rates. Can not reliably determine species. See the top_host_species_additional_stats for a list of probable species from the top returned family " 
    
    for (i in c(1:nrow(hostspeciespercents))) {
      hostspeciespercents$species[i] <- gsub(pattern="s_",x=downscalefam_species$name[i],replacement="")
    }
    
    hostspeciespercents$percentage_reads_total_from_marker <- downscalefam_species$percentage
    
    hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
    
    
    write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
    write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
    
    
    downscalegenus <- "NA"
    downscalegenus2 <- "NA"
    
  } # returns top species from top family and warning indicating multiple species from different genera were potentially the host species
  
  
  if (downscalegenus=="no") {
    
    
    topspeciesgenera1 <- subset(topgenerareturns1,topgenerareturns1$taxRank=="S")
    topspeciesgenera1ordered <- topspeciesgenera1[order(-topspeciesgenera1$cladeReads),]
    
    
    
    topspeciesgenera1top5percent <- subset(topspeciesgenera1ordered,topspeciesgenera1ordered$percentage >=1)
    
    if (nrow(topspeciesgenera1top5percent) >=1) {
    
    hostspecies <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=1))
    hostspeciespercents <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1top5percent),ncol=3))
    colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
    
    
    for (i in c(1:nrow(hostspecies))) {
      hostspecies$V1[i] <- gsub(pattern="s_",x=topspeciesgenera1top5percent$name[i],replacement="")
    }
    
    
    hostspeciespercents$species <- hostspecies$V1
    hostspeciespercents$percentage_reads_total_from_marker <- topspeciesgenera1top5percent$percentage
    
    hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
    
    
    write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
    write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
    }
    
    if (nrow(topspeciesgenera1top5percent) ==0) {
      
      hostspecies <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1ordered),ncol=1))
      hostspeciespercents <- as.data.frame(matrix(nrow=nrow(topspeciesgenera1ordered),ncol=3))
      colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
      
      
      for (i in c(1:nrow(hostspecies))) {
        hostspecies$V1[i] <- gsub(pattern="s_",x=topspeciesgenera1ordered$name[i],replacement="")
      }
      
      
      hostspeciespercents$species <- hostspecies$V1
      hostspeciespercents$percentage_reads_total_from_marker <- topspeciesgenera1ordered$percentage
      
      hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
      
      
      write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
      write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
    }
    
    downscalegenus <- "NA"
    downscalegenus2 <- "NA"
    
  }   # Returns top species from top genus of top family. No warning returned
  
  
}
    
    
    if (downscalefam=="yes") {
      
      
      downscalefam_species <-  subset(topfamreturns,topfamreturns$taxRank=="S")
      
      if(sum(downscalefam_species$percentage >=5) >= 1) {
        downscalefam_speciestop5percent <- subset(downscalefam_species,downscalefam_species$percentage >=5)
      }
      
      else if (sum(downscalefam_species$percentage >=5) == 0) {
        
        downscalefam_speciestop5percent <- downscalefam_species
        
        paste0( " No species found in the top family has over 5% of reads assigned to it. The dataset is likely not resolving enough to determine species for this marker ")
        
      }
      
      hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
      hostspeciespercents <- as.data.frame(matrix(nrow=nrow(downscalefam_speciestop5percent),ncol=3))
      colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
      
      hostspecies[1] <- "Can not reliably determine host species for this marker due to divergence at the family level. See the top_host_species_additional_stats for a list of probable species from the top returned family "
      
      for (i in c(1:nrow(hostspeciespercents))) {
        hostspeciespercents$species[i] <- gsub(pattern="s_",x=downscalefam_speciestop5percent$name[i],replacement="")
      }
      
      hostspeciespercents$percentage_reads_total_from_marker <- downscalefam_speciestop5percent$percentage
      
      hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp <- ((hostspeciespercents$percentage_reads_total_from_marker / sum(hostspeciespercents$percentage_reads_total_from_marker))*100)
      
      
      write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
      write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
      
      
      downscalegenus <- "NA"
      downscalegenus2 <- "NA"
      
      
    } # Returns top species from top family and warning indicating the host family couldn't be reliably returned.
    
  } 
  
  if(nrow(family)==0) {
    
    hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
    hostspeciespercents <- as.data.frame(matrix(nrow=1,ncol=3))
    colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
    
    hostspecies[1] <- "No assignments were relibly inferred for family level or lower " 
    
    hostspeciespercents[1,1] <- NA
    hostspeciespercents[1,2] <- 0
    hostspeciespercents[1,3] <- 0
    
    
  }
  
}


if(nrow(speciesrank)<1) {
  
  hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
  hostspeciespercents <- as.data.frame(matrix(nrow=1,ncol=3))
  colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
  
  hostspecies[1] <- "No assignments were relibly inferred for family level or lower " 
  
  hostspeciespercents[1,1] <- NA
  hostspeciespercents[1,2] <- 0
  hostspeciespercents[1,3] <- 0
  
  write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
  write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
  
  
}
