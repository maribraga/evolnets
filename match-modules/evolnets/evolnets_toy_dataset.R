library(evolnets)
library(bipartite)
library(tidyverse)
library(ape)
library(ggtree)


setwd("~/repos/evolnets/match-modules")

# read trees and character history ----
tree <- read.tree("./evolnets/tree_with_node_labels.tre")
host_tree <- read.tree("./evolnets/host_5tips.tre")

ggt <- ggtree(tree) + geom_tiplab() + geom_nodelab() + theme_tree2() + scale_x_continuous(labels = abs)
revts(ggt)

history <- read_history("./inference/output/out.history.txt")


# get posterior at ages ----
ages <- c(5,2.5,1)
at_ages <- posterior_at_ages(history, ages, tree, host_tree)

samples <- at_ages[[1]]
pps <- at_ages[[2]]

summary_nets_50 <- get_summary_network(pps, 0.5)
summary_nets_50_bin <- get_summary_network(pps, 0.5, weighted = F)


# find modules ----

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

plotModuleWeb(wmod50_5, labsize = 0.6)
plotModuleWeb(wmod50_2.5, labsize = 0.6)
plotModuleWeb(wmod50_1, labsize = 0.6)
