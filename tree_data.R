### Phylogeny analysis

# Read in tree data 
# read.nexus format 
# estimate of branch lengths 
# brain.geomorph.trees.nex

library(ape)
tree <- read.nexus("BrainGeomorph.nex")
plot(tree)

install.packages("tidytree")
library(tidytree)
library(treeplyr)

brain.data <- make.treedata(tree, S1.GPA)
