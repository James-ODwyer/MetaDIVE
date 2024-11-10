library(argparse)

parser <- ArgumentParser(description= 'Summarising host id files')

parser$add_argument('--inputCO1', '-a', help= 'Input taxonomy file')
parser$add_argument('--inputLSU', '-b', help= 'Input taxonomy file')
parser$add_argument('--inputSSU', '-c', help= 'Input taxonomy file')
parser$add_argument('--programdir', '-d', help= 'Base directory')
parser$add_argument('--name', '-e', help= 'Base directory')
parser$add_argument('--outdir', '-f', help= 'output directory')
parser$add_argument('--CO1weighting', '-x', help= 'Weighting value for CO1')
parser$add_argument('--LSUweighting', '-y', help= 'Weighting value for LSU')
parser$add_argument('--SSUweighting', '-z', help= 'Weighting value for SSU')

xargs<- parser$parse_args()

inputCO1 <- read.table(xargs$inputCO1,sep="\t", header=TRUE,row.names = NULL)

inputSSU <- read.table(xargs$inputSSU,sep="\t", header=TRUE,row.names = NULL)

inputLSU <- read.table(xargs$inputLSU,sep="\t", header=TRUE,row.names = NULL)

CO1weight <- xargs$CO1weighting
CO1weight <- as.numeric(CO1weight)
LSUweight <- xargs$LSUweighting
LSUweight <- as.numeric(LSUweight)
SSUweight <- xargs$SSUweighting
SSUweight <- as.numeric(SSUweight)

NAMES <- xargs$name

basepath <- xargs$programdir
resultspath <- xargs$outdir

outtablespath <- paste0(basepath,resultspath)


inputCO1$value <- (CO1weight * inputCO1$percentage_reads_total_from_marker)
inputSSU$value <- (SSUweight * inputSSU$percentage_reads_total_from_marker)
inputLSU$value <- (LSUweight * inputLSU$percentage_reads_total_from_marker)

allinputs <- rbind(inputCO1,inputSSU,inputLSU)



summaryvalues <- aggregate(allinputs$value, by=list(Species=allinputs$species), FUN=sum)

if(nrow(summaryvalues) >=1) {
summaryvaluesordered <- summaryvalues[order(-summaryvalues$x),]


top_host_species <- summaryvaluesordered$Species[1]

write.table(top_host_species,file=(paste0(outtablespath,NAMES,"_top_host_species_overall.txt")), sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(summaryvaluesordered,file=(paste0(outtablespath,NAMES,"_top_host_scaled_host_info.txt")),sep="\t",quote = FALSE, row.names = FALSE, col.names = TRUE)

}

if(nrow(summaryvalues) ==0) {
top_host_species <- "No species assigned due to insufficient reads in any marker"
write.table(top_host_species,file=(paste0(outtablespath,NAMES,"_top_host_species_overall.txt")), sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(summaryvalues,file=(paste0(outtablespath,NAMES,"_top_host_scaled_host_info.txt")),sep="\t",quote = FALSE, row.names = FALSE, col.names = TRUE)


}