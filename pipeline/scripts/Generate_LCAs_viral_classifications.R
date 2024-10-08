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

xargs<- parser$parse_args()

inputtaxa <- xargs$inputtax

NAMES <- xargs$name

basepath <- xargs$programdir
resultspath <- xargs$outputpath

outtablespath <- paste0(basepath,resultspath)

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
      dragY = TRUE,
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 12,
      linkGradient = TRUE,
      nodeShadow = TRUE,
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

if (!(is.null(myfile))) {

sankeynetwork <- build_sankey_network(myfile,taxRanks = c("D","K","F","G","S"),zoom = FALSE)

paste0("saving output figure to the following ", outtablespath,NAMES,"_viral_taxonomic_classifications_of_raw_reads", ".html")

sankeynetwork %>% sankeyD3::saveNetwork(file =paste0(outtablespath,NAMES,"_viral_taxonomic_classifications_of_raw_reads", ".html"))

}


if ((is.null(myfile))) {
  
  paste0("Kraken returned zero assignments to species for any virus. Report cannot be generated")
}

# Had the host id part in the original barcode script here. 
# Not sure if there is merit in creating an identify genus/families undefined here 
# but that would be a good spot for it here as it is largely similar in how I would identify them