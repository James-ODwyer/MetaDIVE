library(argparse)

parser <- ArgumentParser(description= 'Summarising results filtering and assembly')

parser$add_argument('--genesannot', '-a', help= 'the annotated genes file created by Dram and formatted in python')
parser$add_argument('--outfile', '-o', help= 'Output file and path for contig ids and target bp regions')


xargs<- parser$parse_args()

resultsfile <- xargs$outfile

genes_annot <- readLines(xargs$genesannot)


grep(pattern = ".*Dbxref.*viral",x = genes_annot) -> rowidx


genes_annot_viral <- genes_annot[rowidx]




genes_annot_viralstrings <- strsplit(genes_annot_viral,split = "\t")

genesannot_viral_table <- as.data.frame(matrix(nrow=length(genes_annot_viral),ncol=3))


for ( i in c(1:length(genes_annot_viral))) {
  
  
  genes_annot_viralstrings[[i]][1] -> genesannot_viral_table[i,1]
  genes_annot_viralstrings[[i]][4] -> genesannot_viral_table[i,2]
  genes_annot_viralstrings[[i]][5] -> genesannot_viral_table[i,3]
  
}

genesannot_viral_table_uniq <- unique(genesannot_viral_table)


genesannot_viral_table_uniq$V1 <- gsub(pattern = "__.*",replacement = "",genesannot_viral_table_uniq$V1)


write.table(genesannot_viral_table_uniq,file=(paste0(resultsfile)),sep="\t",row.names=FALSE, quote=FALSE,col.names = FALSE)