library(ggtree)

tree  <- rcoal(6)
tree$tip.label <- paste0("T",1:6)
tree$edge.length
tree.height <- node.depth.edgelength(tree)[1]
tree$edge.length <- tree$edge.length*10/tree.height
ggtree(tree, ladderize = F) + geom_tiplab() +theme_tree2()

write.tree(tree, "tree_6tips.tre")

host <- rcoal(5)
host$tip.label <- paste0("H",1:5)
plot(host)
write.tree(host, "host_5tips.tre")

tnodes <- ape::read.tree("./inference/tree_with_node_labels.tre")


