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

parser <- ArgumentParser(description= 'Summarising results raws and contigs combined')

parser$add_argument('--inputbasedir', '-a', help= 'program run directory')
parser$add_argument('--input_resultscsvsmain ', '-b', help= 'master directory of all sample folders which contain output csvs')
parser$add_argument('--outputpath', '-c', help= 'summary_results2 output directory')
parser$add_argument('--minreadthreshold', '-d', help= 'Minimum reads required to report in summary statistics')

xargs<- parser$parse_args()

# define base parameters and parameter variables for which analyses were run
basepath <- xargs$inputbasedir
resultspath <- xargs$outputpath
inputdir <- xargs$input_resultscsvsmain
outtablespath <- paste0(basepath,resultspath)
intablespath <- paste0(basepath,inputdir)



readsfilterthreshold <- xargs$minreadthreshold


  summary_returned_csv <- list.files(path = intablespath, pattern = ".*/*_virusall_sums.csv", all.files = FALSE,
                                               full.names = TRUE, recursive = TRUE,
                                               ignore.case = FALSE, include.dirs = FALSE)

all_data <- list()

# Loop through each file path in summary_returned_csv
for (file in summary_returned_csv) {
  
  # Check the number of rows in the file (count lines)
  num_lines <- length(readLines(file))
  
  # Only read the file if it has more than one line (header + data)
  if (num_lines > 1) {
    
    # Read each CSV file
    data <- read.csv(file)
    
    # Remove the column "max_percent_ident.1" if it exists
    if ("max_percent_ident.1" %in% names(data)) {
      data <- data[ , !names(data) %in% "max_percent_ident.1"]
    }
    
    # Extract the sample name by splitting the path and getting the second-to-last component
    sample_name <- basename(dirname(file))
    
    # Add the sample name as a new column
    data$Sample <- sample_name
    
    # Append the dataframe to the list
    all_data[[file]] <- data
  }
}

# Combine all dataframes into one
viruses_all_reads <- do.call(rbind, all_data)

# remove really long rownames as they make it impossible to look at 
row.names(viruses_all_reads) <- NULL

#subset to the minimum viral count threshold
viruses_reads_threshold <- subset(viruses_all_reads, viruses_all_reads$total_reads_assigned >=readsfilterthreshold)

viruses_reads_threshold <- viruses_reads_threshold[order(-viruses_reads_threshold$total_reads_assigned),]

if (nrow(viruses_reads_threshold)>=10) {
  viruses_top10 <- viruses_reads_threshold[1:10,]
}
if (nrow(viruses_reads_threshold)<10) {
  viruses_top10 <- viruses_reads_threshold
}



viruses_top10$total_reads_assigned <- as.numeric(as.character(viruses_top10$total_reads_assigned))


plot <- ggplot(viruses_top10, aes(x=Sample, y=total_reads_assigned, fill=subspecies)) +
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


ggsave(paste0(outtablespath,"Top_viral_hits_combined_raws_plus_contigs.pdf"), dpi = 330)

plot 

ggsave(paste0(outtablespath,"Top_viral_hits_combined_raws_plus_contigs.png"), dpi = 330)




write.table(viruses_reads_threshold,file=(paste0(outtablespath,"Top_viral_hits_combined_raws_plus_contigs.txt")),sep="\t",row.names=FALSE,col.names=colnames(viruses_reads_threshold))



