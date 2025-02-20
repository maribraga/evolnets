% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index_at_ages.R
\name{index_at_ages_summary}
\alias{index_at_ages_summary}
\title{Network structure indices across ancestral summary networks and extant network}
\usage{
index_at_ages_summary(
  summary_networks,
  index,
  ages = NULL,
  nnull = 100,
  use_future = FALSE
)
}
\arguments{
\item{summary_networks}{List of ancestral networks at given ages and extant network. Usually from
the output of \code{get_summary_networks}, see example.}

\item{index}{Index to be calculated for each ancestral network. "NODF" to calculate nestedness or
"Q" to calculate modularity. For weighted networks, "weighted NODF" will be calculated.}

\item{ages}{Vector of ages (time points in the past) of ancestral networks. By default, uses all
ages present in \code{summary_networks}.}

\item{nnull}{Number of null networks to generate to calculate the z-score. Default is 100.}

\item{use_future}{Parallelize with package \code{future}? Logical. Only applicable to modularity calculation.}
}
\value{
A data.frame of z-scores and p-values across networks.
}
\description{
Calculate z-scores for nestedness (NODF) or modularity (Q) for each ancestral network at time points in
the past and the extant network, based on null networks where all interactions have the same probability. By calculating
z-scores, we can compare ancestral networks at different ages.
}
\examples{
# read data that comes with the package
data_path <- system.file("extdata", package = "evolnets")
tree <- read_tree_from_revbayes(paste0(data_path,"/tree_pieridae.tre"))
host_tree <- ape::read.tree(paste0(data_path,"/host_tree_pieridae.phy"))
history <- read_history(paste0(data_path,"/history_thin_pieridae.txt"))

# get ancestral networks at ages in the past
ages <- c(60, 50, 40, 0)
at_ages <- posterior_at_ages(history, ages, tree, host_tree)
summary_networks <- get_summary_networks(at_ages, threshold = 0.5, weighted = TRUE)

# calculate nestedness of ancestral and extant networks
Nz <- index_at_ages_summary(summary_networks, index = "NODF")

# calculate modularity of ancestral and extant networks with parallelization (slower)
# Qz <- index_at_ages_summary(summary_networks, index = "Q", use_future = TRUE)
}
