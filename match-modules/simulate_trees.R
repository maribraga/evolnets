library(ggtree)
library(ape)

tree  <- rcoal(20)
tree$tip.label <- paste0("T",1:Ntip(tree))
tree$edge.length
tree.height <- node.depth.edgelength(tree)[1]
tree$edge.length <- tree$edge.length*100/tree.height
ggtree(tree, ladderize = F) + geom_tiplab() +theme_tree2()

write.tree(tree, paste0("./inference/data/tree_",Ntip(tree),"tips.tre"))

host <- rcoal(10)
host$tip.label <- paste0("H",1:Ntip(host))
plot(host)
write.tree(host, paste0("./inference/data/host_",Ntip(host),"tips.tre"))

# data for inference
matrix <- read.csv("toy_example.csv", header = T, row.names = 1)
write.nexus.data(matrix, "./inference/data/toy_data.nex")

