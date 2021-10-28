library(evolnets)
library(bipartite)
library(tidyr)
library(dplyr)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)
library(gridExtra)
library(gridBase)
library(grid)
library(restorepoint)

setwd("~/repos/evolnets/match-modules")
#setwd("~/projects/evolnets/match-modules")

# read trees and character history ----

# 6 symbionts and 5 hosts
#host_tree <- read.tree("./evolnets/host_5tips.tre")
#tree <- read.tree("./evolnets/.tre")  # find the right tree


# 20 symbionts and 10 hosts
host_tree <- read.tree("./evolnets/host_10tips.tre")
treeRev <- read.beast.newick("./evolnets/tree20_Rev.tre")
net <- as.matrix(read.csv("./evolnets/incidence_matrix_20_10.csv", header = T, row.names = 1))
tree <- treeRev@phylo

tree$node.number <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))

index_node <- treeRev@data %>%
  mutate(node = as.numeric(node)) %>%
  arrange(node)

indices <- index_node %>% filter(node > Ntip(tree)) %>% pull(index)
names(indices) <- NULL

tree$node.label <- paste0("Index_",indices)

plot(tree, show.node.label = T)

ggt <- ggtree(tree, ladderize = F) +
  geom_tiplab(size = 3) +
  geom_nodelab(size = 3) +
  theme_tree2() +
  scale_x_continuous(labels = abs)
ggt <- revts(ggt)
ggt + geom_vline(xintercept = c(-50,-40,-30,-20,-10), lty=2, col="red")
ggt

# get posterior at ages ----
history <- read_history("./inference/output/out.history.txt")

ages <- c(50,40,30,20,10,0)
at_ages <- posterior_at_ages(history, ages, tree, host_tree)

samples <- at_ages[[1]]
pps <- at_ages[[2]]

summary_networks <- get_summary_network(pps, ages, 0.5)
#summary_nets_50_bin <- get_summary_network(pps, ages, 0.5, weighted = F)


# find modules ----

# at once
all_mod <- modules_across_ages(summary_networks, tree)
matched_mod <- all_mod[[2]]
#write.csv(matched_mod,"./evolnets/match_mod_output.csv", row.names = F)
unmatched_modwebs <- all_mod[[1]][[2]]
unmatched_modules <- all_mod[[1]][[1]]

par(mfrow=c(2,3))
for (i in 1:length(ages)) {
    plotModuleWeb(unmatched_mod[[i]], labsize=0.6)
}

# plot asr
at_nodes <- posterior_at_nodes(history, nodes = (Ntip(tree)+1):(Nnode(tree)+Ntip(tree)), host_tree)
modules <- all_mod %>%
  filter(age == 0) %>%
  mutate(module = as.character(match(module_name, LETTERS)),
         type = case_when(grepl("T", name) ~ "symbiont",
                          grepl("H", name) ~ "host")) %>%
  select(name,module,type)

p <- plot_ancestral_states(tree,at_nodes,modules,threshold = 0.5)
p + geom_vline(xintercept = c(-50,-40,-30,-20,-10), lty=2, col="red")

plot_module_matrix2(net,at_nodes, tree,host_tree,modules,threshold = 0.7)

# to plot unmatched modules
# list_all_mod <- modules_from_summary_networks(summary_networks)
#
# par(mfrow = c(2,3))
# for(p in 1:length(ages)){
#   plotModuleWeb(list_all_mod[[2]][[p]], labsize=0.6)
# }
#
# unmatched_modules <- list_all_mod[[1]]
#saveRDS(unmatched_modules, "./evolnets/unmatched_modules_20_10.rds")
unmatched_modules <- readRDS("./evolnets/unmatched_modules_20_10.rds")

# match modules
all_mod <- match_modules(summary_networks, unmatched_modules, tree) #, module_strength = "mean")

restore.point("point1", to.global = F)






## ---
all_wmod50 <- tibble()

for(i in 1:length(summary_nets_50)){
  set.seed(5)
  wmod <- computeModules(summary_nets_50[[i]])
  assign(paste0("wmod50_",ages[i]),wmod)
  wmod_list <- listModuleInformation(wmod)[[2]]
  nwmod <- length(wmod_list)

  for(m in 1:nwmod){
    members <- unlist(wmod_list[[m]])
    mtbl <- tibble(name = members,
                   age = rep(ages[i], length(members)),
                   original_module = rep(m, length(members)))

    all_wmod50 <- bind_rows(all_wmod50, mtbl)
  }
}

# plot trees and modules
mod_list = list(  wmod50_5, wmod50_2.5, wmod50_1, wmod50_0)
tree_list = list()
for (i in 1:length(ages)) {
    tree_list[[i]] = ggt + geom_vline(xintercept = -ages[i], lty=2, col="red")
}

pdf( "./evolnets/toy_modules.pdf", height=5, width=12 )
plot.new()
#grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 4)))

for (i in 1:length(ages)) {
    # ggplot
    pushViewport(viewport(layout.pos.col = i, layout.pos.row=1))
    print(tree_list[i], newpage = FALSE)
    popViewport()

    #Draw bsae plot
    pushViewport(viewport(layout.pos.col = i, layout.pos.row=2))
    par(fig = gridFIG(), new = TRUE)
    plotModuleWeb(mod_list[[i]], labsize=0.6)
    popViewport()
}
dev.off()


### Function: for a given node/branch (x) at a given time (t_{i}),
#             find all daughter nodes (y) at the next time step (t_{i+1}),
#             and identify all modules for those nodes at that time;
#             then have parent nodes "reverse-inherit" the modules of
#             the child nodes
#
#             if an ancestral node/branch has relationships
#             with two daughter modules, give the daughter module name to the
#             ancestor with the highest total probability of interactions;
#             give the other ancestors(s) new module names
#
#             if daughter modules merge ancestrally, then retain the name
#             for the largest extant module


test <- match_modules(all_wmod50, ages, tree)


# get all nodes between 2 time slices

ages <- sort(ages)
node.depth.edgelength(tree)


