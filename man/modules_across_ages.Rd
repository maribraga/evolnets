% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules.R
\name{modules_across_ages}
\alias{modules_across_ages}
\title{Evolution of modules: find and match modules across time slices}
\usage{
modules_across_ages(summary_networks, tree, extant_modules = NULL)
}
\arguments{
\item{summary_networks}{List of reconstructed summary networks for each age (output from \code{get_summary_network()}).}

\item{tree}{The phylogeny of the symbiont clade (e.g. parasite, herbivore), a \code{phylo} object.}

\item{extant_modules}{A \code{moduleWeb} object defining the modules in the extant network.}
}
\value{
A list with:
\enumerate{
\item A list of 2 elements: 1.1) a data frame containing the module information for each node at each
network, 1.2) a data frame of correspondence between the original and the matched module names for each network;
\item A list of 2 elements: 2.1) a data frame containing the module membership of each node at each age before matching, 2.2) a list of \code{moduleWeb} objects for each age.
}
}
\description{
Evolution of modules: find and match modules across time slices
}
\examples{
# read data that comes with the package
data_path <- system.file("extdata", package = "evolnets")
tree <- read_tree_from_revbayes(paste0(data_path,"/tree_pieridae.tre"))
host_tree <- ape::read.tree(paste0(data_path,"/host_tree_pieridae.phy"))
history <- read_history(paste0(data_path,"/history_thin_pieridae.txt"), burnin = 0)

# get ancestral summary networks
ages <- c(60, 50, 40, 0)
at_ages <- posterior_at_ages(history, ages, tree, host_tree)
summary_networks <- get_summary_networks(at_ages, threshold = 0.5, weighted = TRUE)

# find and match modules across ancestral and extant networks
all_mod <- modules_across_ages(summary_networks, tree)
}
