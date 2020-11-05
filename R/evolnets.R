#' evolnets: a package for summarizing inference of network evolution
#'
#' The evolnets package provides three categories of important functions:
#' rates, ancestral states and samples.
#'
#' @section Rates:
#' These functions are used to calculate effective rates of host-repertoire
#'     evolution: effective_rate( ), count_events( ), rate_gl( ), count_gl( ).
#'
#' @section Ancestral states:
#' These functions are used to calculate the posterior probabilities of
#'     host-parasite interactions at internal nodes of the parasite tree or at
#'     specific time points in the past: posterior_at_nodes( ), posterior_at_ages( ).
#'
#' @section Samples:
#' These functions perform calculations for each sampled host-parasite
#'     network during MCMC: samples_at_ages( ), Q_posterior_at_ages( ),
#'     NODF_posterior_at_ages( ).
#'
#' @section More information:
#' An example of how to use the package can be found in the vignette [**Introduction to evolnets**](evolnets.html).
#'
#'
#' @docType package
#' @name evolnets
NULL
#> NULL
