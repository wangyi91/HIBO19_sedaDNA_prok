library(ape) # Load the ape package
library(picante)
library(stringr)
library(readr)
library(dplyr)
library(treeio)

numbers_only <- function(x) {suppressWarnings(!is.na(as.numeric(x)))}

# Data preparation
# 1) taxa list and Nread output from julia. Path: 
ftax="Bacteria";
ftax="Archaea";
# 2) .newick trees from PhyloT. Path:
tree.dir <- paste("./data/", ftax, sep='')

# process tree per sample
for (infile in dir(tree.dir, pattern="_HB_[0-9][0-9].newick$")) {
  lib = str_extract(infile, "HB_[0-9]+")
  tree <- read.tree(paste(tree.dir, "/", infile, sep=""))
 
  # Clean tree by rm bootstrap values
  tree$node.label <- ifelse(numbers_only(tree$node.label), "", tree$node.label) # Erase the bootstrap values from the phylo object
  write.tree(tree, file = paste(tree.dir,'/', gsub("\\.newick","",infile), ".modified.newick", sep="")) # Save it
  
}

# Now send trees back to calculate PDs using ./MicrobeProfiling/phylo_diversity.jl.


