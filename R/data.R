#' Data for Pieridae butterflies and their host plants from Braga et al. 2020
#'
#' @format `history`: a data frame with 58,891 rows and 14 columns:
#' \describe{
#'   \item{iteration}{MCMC iteration}
#'   \item{posterior, likelihood, prior}{Probabilities}
#'   \item{node_index}{Index of the tree node}
#'   \item{branch_start_time, branch_end_time}{Time when the branch starts and ends}
#'   \item{start state, end state}{Host repertoire at the start and end of the given branch}
#'   \item{transition_time}{Time when event (host gain or loss) happened}
#'   \item{parent_index, child1_index, child2_index}{Index of parent and child nodes}
#' }
#' @details The data frame contains MCMC samples from the analysis of host-repertoire
#'  evolution across Pieridae, done in RevBayes
#'
#' @source Braga et al. 2020 BioRxiv
"history"

#' Host plants tree
#' @format `host_tree`: a phylogenetic tree of angiosperm families with 50 tips.
#' @source Braga et al. 2020 BioRxiv
"host_tree"

#' Pieridae tree
#' @format `tree`: a phylogenetic tree of Pieridae genera with 66 tips.
#' @source Braga et al. 2020 BioRxiv
"tree"


