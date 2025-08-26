

if (!require("d3treeR", character.only = TRUE)) {
  devtools::install_github("timelyportfolio/d3treeR",upgrade="never",dependencies=FALSE)
  library("d3treeR")
} else {
  library("d3treeR")
}


if (!require("d3treeR", character.only = TRUE)) {
  devtools::install_github("timelyportfolio/d3treeR",upgrade="never")
  library("d3treeR")
} else {
  library("d3treeR")
}