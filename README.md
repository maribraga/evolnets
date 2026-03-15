
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evolnets

<!-- badges: start -->

<!-- badges: end -->

RevBayes and TreePPL offer models to infer host repertoire evolution,
but no tools to parse the outputs. *evolnets* has the necessary tools to
reconstruct ancestral ecological networks based on posterior
probabilities of interactions.

## Installation

You can install evolnets like so:

``` r
if(!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
  library(devtools)
} else {
 library(devtools)
}

devtools::install_github("evonetslab/evolnets")
```

## About

The evolnets package provides three categories of important functions:
rates, ancestral states and samples.

- Rates: these functions are used to calculate effective rates of
  host-repertoire evolution: `effective_rate( )`, `count_events( )`,
  `rate_gl( )`, `count_gl( )`.

- Ancestral states: these functions are used to calculate the posterior
  probabilities of host-parasite interactions at internal nodes of the
  parasite tree or at specific time points in the past:
  `posterior_at_nodes( )`, `posterior_at_ages( )`.

- Samples: these functions perform calculations for each sampled
  host-parasite network during MCMC: `samples_at_ages( )`,
  `Q_posterior_at_ages( )`, `NODF_posterior_at_ages( )`.

See the full documentation at the [evolnets’
website](https://evonetslab.github.io/evolnets/).
