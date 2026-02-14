library(ape)
library(picante)
library(stringr)
library(readr)
library(dplyr)
library(treeio)

numbers_only <- function(x) {suppressWarnings(!is.na(as.numeric(x)))}

ftax="Bacteria"; 
ftax="Archaea";
# read .newick trees from PhyloT. Path:
tree.dir <- paste("./data/", ftax, sep='')

# read the tree using all taxa to get a list to ranknames
tree <- read.tree(paste(tree.dir, "/","comm_",ftax,".newick", sep="")) #

# Clean tree by rm bootstrap values
tree$node.label <- ifelse(numbers_only(tree$node.label), "", tree$node.label) # Erase the bootstrap values from the phylo object

# save tree
write.tree(tree, file = paste(tree.dir, "/","comm_",ftax,".modified.newick", sep="")) # Save it



# Now send trees back to julia





