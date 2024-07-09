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
  
  hostspecies <- as.data.frame(matrix(nrow=1,ncol=1))
  hostspeciespercents <- as.data.frame(matrix(nrow=1,ncol=3))
  colnames(hostspeciespercents) <- c("species", "percentage_reads_total_from_marker","percentage_reads_scaled_to_all_returned_top_sp")
  
  hostspeciespercents$species <- "NA"
  hostspeciespercents$percentage_reads_total_from_marker <- 0
  hostspeciespercents$percentage_reads_scaled_to_all_returned_top_sp  <- 0
  
  write.table(hostspecies,file=(paste0(outtablespath,NAMES,"_top_host_species_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE, col.names = FALSE)
  write.table(hostspeciespercents,file=(paste0(outtablespath,NAMES,"_top_host_species_additional_stats_",marker, ".txt")),sep="\t", quote =FALSE, row.names=FALSE,col.names = TRUE)
