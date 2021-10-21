library(ggtree)
library(ape)

ntips <- 20
tree  <- rcoal(ntips)
tree$tip.label <- paste0("T",1:ntips)
tree$edge.length
tree.height <- node.depth.edgelength(tree)[1]
tree$edge.length <- tree$edge.length*100/tree.height
ggtree(tree, ladderize = F) + geom_tiplab() +theme_tree2()

write.tree(tree, paste0("./inference/data/tree_",Ntip(tree),"tips.tre"))

nhosts <- 44
host <- rcoal(nhosts)
host$tip.label <- paste0("H",1:nhosts)
plot(host)
write.tree(host, paste0("./inference/data/host_",Ntip(host),"tips.tre"))

# data for inference
matrix <- read.csv("toy_example.csv", header = T, row.names = 1)
write.nexus.data(matrix, "./inference/data/toy_data.nex")

