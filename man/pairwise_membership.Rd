% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules.R
\name{pairwise_membership}
\alias{pairwise_membership}
\title{Calculate frequency that pairs of nodes fall within the same module}
\usage{
pairwise_membership(mod_samples, ages, edge_list = TRUE)
}
\arguments{
\item{mod_samples}{Output from \code{modules_from_samples()}}

\item{ages}{Vector of network ages}

\item{edge_list}{Should output be an edge list?}
}
\value{
An edge list with the frequency that each pair of nodes in the network is placed in the same module across network samples.
}
\description{
Calculate frequency that pairs of nodes fall within the same module
}
\examples{
\dontrun{
pairwise_membership(mod_samples, c(0))
}
}
