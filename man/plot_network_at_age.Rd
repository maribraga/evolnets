% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_network_at_age}
\alias{plot_network_at_age}
\title{Plot one ancestral network with module information at a given time}
\usage{
plot_network_at_age(
  subtree,
  tip_data,
  tgraph,
  module_levels,
  colors = NULL,
  tree,
  age,
  weighted = TRUE,
  weight_range = c(0, 1),
  two_state = FALSE,
  node_size = NULL,
  level_names = c("Host", "Symbiont")
)
}
\arguments{
\item{subtree}{a \code{phylo} object of the original tree sliced at a given time in the past.}

\item{tip_data}{a \code{data.frame} containing the module information for each tip in the subtree.}

\item{tgraph}{a \code{tbl_graph} containing the nodes and edges of the ancestral network and the
module information for each node.}

\item{module_levels}{Order in which the modules should be organized. Affects which color each
module will be assigned.}

\item{colors}{Color vector used to plot module information.}

\item{tree}{The phylogeny of the symbiont clade (e.g. parasites, herbivores). Object of class
\code{phylo}.}

\item{age}{Age of the ancestral network to be plotted as the tittle.}

\item{weighted}{Whether the network should have weighted edges.}

\item{weight_range}{The range of weights conscidered for the width of the edges.}

\item{two_state}{Whether the width of the edges should reflect the state (instead of the
posterior probability).}

\item{node_size}{Size of the nodes in every network. If NULL, the default size is plotted.}

\item{level_names}{Optional string vector with two elements used by legend: name of the clade
modeled as "host" and the name of the clade modeled as "symbiont".}
}
\value{
An assembly of plots, of class \code{patchwork}.
}
\description{
Plot one ancestral network with module information at a given time
}
\examples{
\dontrun{
plot_network_at_age(subtree, tip_data, tgraph, module_levels, colors, tree, age)
}
}
