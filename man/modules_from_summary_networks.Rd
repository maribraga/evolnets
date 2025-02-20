% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules.R
\name{modules_from_summary_networks}
\alias{modules_from_summary_networks}
\title{Identify modules for each summary network at each age}
\usage{
modules_from_summary_networks(summary_networks, extant_modules = NULL)
}
\arguments{
\item{summary_networks}{List of reconstructed summary networks for each age (output from \code{get_summary_network()}).}

\item{extant_modules}{A \code{moduleWeb} object defining the modules in the extant network.}
}
\value{
A list of 2 elements: 1) a data frame containing the module membership
of each node at each age before matching; 2) a list of \code{moduleWeb} objects for each age.
}
\description{
This function is called within \code{modules_across_ages()}.
}
\examples{
\dontrun{
 unmatched_modules <- modules_from_summary_networks(summary_networks)
}
}
