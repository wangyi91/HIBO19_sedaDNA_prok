library(ape) # Load the ape package
library(picante)
library(stringr)
library(readr)
library(dplyr)
library(treeio)

numbers_only <- function(x) {suppressWarnings(!is.na(as.numeric(x)))}

# Data preparation
# 1) taxa list and Nread output from julia. Path: 
file.dir <- "/Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/MicrobeProfiling/output"

ftax="Bacteria";
ftax="Archaea";
# 2) .newick trees from PhyloT. Path:
tree.dir <- paste("./data/", ftax, sep='')

# read the tree using all taxa to get a list to ranknames
alltree <- read.tree(paste(tree.dir, "/", ftax, add3, "_all_50.newick", sep=""))
ranknames = unique(alltree$node.label) %>% as.data.frame %>% filter(startsWith(.,"p__"))
write_csv(ranknames, paste(tree.dir,"/ranknames_",ftax, add3,".csv",sep=''), quote="none")


for (infile in dir(tree.dir, pattern="_HB_[0-9][0-9].newick$")) {
  lib = str_extract(infile, "HB_[0-9]+")
  tree <- read.tree(paste(tree.dir, "/", infile, sep=""))
 
  # Clean tree by rm bootstrap values
  tree$node.label <- ifelse(numbers_only(tree$node.label), "", tree$node.label) # Erase the bootstrap values from the phylo object
  #tree <- collapse.singles(tree)
  write.tree(tree, file = paste(tree.dir,'/', gsub("\\.newick","",infile), ".modified.newick", sep="")) # Save it
  
  # subset tree at phylum and save subtree
  #for (rank in ranknames$.) {
   # tryCatch({
    #subtree <- tree_subset(tree, rank, levels_back = 0)
    #write.tree(subtree, file = paste(tree.dir,'/', gsub("\\.newick","",infile), ".", rank,".subtree.newick", sep=""))
    #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # ignore if a sample does not have this phylum of taxa
  #}
}




# Now send trees back to julia to calc PDs.






# use picante to calc PD
#for(infile in dir(file.dir, pattern="\\.newick$")) {
  lib = str_extract(infile, "HB_[0-9]+")
  
  tree <- read.tree(paste(tree.dir, "/", infile, sep=""))
  ab <- read_delim(paste(file.dir, "/tax_Nread_filtered_", lib, "_100.csv", sep=""), delim =',', col_names =c("tax","N_reads"))
  ab$tax = gsub(" ","_", ab$tax)
  
  ab <- ab %>% as.data.frame %>% `row.names<-`(.$tax) %>% select(-tax) %>% t
  
  # Clean tree by rm internal nodes and bootstrap values
  tree$node.label <- ifelse(numbers_only(tree$node.label), "", tree$node.label) # Erase the bootstrap values from the phylo object
  #tree <- collapse.singles(tree)
  #write.tree(tree, file = paste(tree.dir, gsub("\\.newick","",infile), ".modified.newick", sep="")) # Save it
  
  comm <- match.phylo.comm(phy = tree, 
                           comm = ab)$comm
  phy <- match.phylo.comm(phy = tree, 
                          comm = ab)$phy
  
  pd(ab, phy, include.root=TRUE)
}



