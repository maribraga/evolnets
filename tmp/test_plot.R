library(tidyverse)
library(ape)
library(ggtree)
library(treeio)

# extant network
net <- read.csv("~/repos/pieridae_hostrep/data/incidence_pieridae.csv", header = TRUE, row.names = 1) %>%
  as.matrix()

# trees
data(tree) # this one already has correct node labels
data(host_tree)

# modules
modules_across_ages <- read.csv("~/repos/pieridae_hostrep/networks/all_wmod50_bl1.csv", header = TRUE, stringsAsFactors = F)
modules <- modules_across_ages %>%
  filter(age == 0) %>%
  select(name, module) %>%
  mutate(module = sub('.', '', module),                            # remove 'M' from each module
         type = case_when(name %in% tree$tip.label ~ "symbiont",
                          name %in% host_tree$tip.label ~ "host"))

# samples_at_nodes
data(history)

nodes <- 67:131
samples_at_nodes <- posterior_at_nodes(history, nodes, host_tree)

# plot everything
plot_module_matrix2(net, samples_at_nodes, tree, host_tree)

plot_module_matrix2(net, samples_at_nodes, tree, host_tree, modules)
