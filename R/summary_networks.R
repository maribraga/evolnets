#' Get summary networks from posterior probabilities of interactions of extant lineages at given
#' time points in the past
#'
#' @param at_ages List of lists, with samples and posterior probabilities for each interaction of
#'   extant lineages at given ages. Usually calculated `posterior_at_ages`, see example.
#' @param threshold Probability threshold to include an interaction in the network. Interactions
#'   with posterior `probability < threshold` will be dropped.
#' @param ages Vector of ages (time points in the past) at which samples were retrieved. By default,
#'   uses all ages present in `pp_at_ages`.
#' @param weighted Logical. Use posterior probabilities as interaction weights?
#' @param type One of `'states'` or `'repertoires'`. If `'states'`, will plot the presence of a
#'   state when its posterior probablity is higher than `threshold`. If `'repertoires'`, will plot
#'   the same but for the given `repertoire`.
#' @param state Which state? Default is 2. For analyses using the 3-state model, choose `1`, `2` or
#'   both using `c(1, 2)` (where 1 is a potential host and 2 an actual host). Only used if `type` is
#'   `'states'`.
#' @param repertoire Either the `'realized'` repertoire which is defined as state 2, or the
#'   `'fundamental'` repertoire (default) which is defined as having any state (usually 1 or 2), and
#'   in the 3-state model includes both actual and potential hosts.
#'
#' @return A list of incidence matrices (summary network) for each time slice in `ages`.
#' @export
#'
#' @examples
#' data(tree)
#' data(host_tree)
#' data(history)
#'
#' ages <- c(60, 50, 40, 0)
#' at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#' weighted_net_50 <- get_summary_network(at_ages, threshold = 0.5, weighted = TRUE)
#' binary_net_90 <- get_summary_network(at_ages, threshold = 0.9, weighted = FALSE)
get_summary_network <- function(
    at_ages, threshold, ages = NULL, weighted = TRUE, type = "states", state = 2,
    repertoire = 'fundamental'
){

  if (!is.list(at_ages)) stop('`pp_at_ages` should be a list.')
  if (!is.numeric(threshold) || !(threshold > 0 & threshold <= 1)) {
    stop('`pt` should be a numeric value between 0 and 1.')
  }
  if (!is.null(ages) && !is.numeric(ages)) stop('`ages` should be a numeric vector or NULL.')
  if (!is.logical(weighted)) stop('`logical` should be a logical value (TRUE/FALSE).')
  if (!(type %in% c('states', 'repertoires') & length(type) == 1)) {
    stop("`type` should be either 'states' or 'repertoires'.")
  }
  if (!all(as.character(state) %in% dimnames(at_ages[['post_states']][[1]])[[3]])) {
    stop("`state` should be a vector and have valid states occuring in `at_ages`")
  }
  if (!(repertoire %in% c('fundamental', 'realized') & length(repertoire) == 1)) {
    stop("`repertoire` should be either 'fundamental' or 'realized'.")
  }
  if (threshold <= 0.5 & type == 'states' & length(state) > 1) {
    stop('When analyzing multiple states, threshold should be > 0.5.')
  }

  # pick the right states or repertoires
  if (type == 'states') {
    if (weighted & (length(state) > 1)) stop('Multiple states are not supported for binary networks.')
    age_list <- at_ages$post_states
    age_list <- lapply(age_list, function(a) a[, , as.character(state)])
  } else {
    age_list <- at_ages$post_repertoires
    age_list <- lapply(age_list, function(a) a[, , as.character(repertoire)])
  }

  if (!is.null(ages) && (length(at_ages) != length(age_list))) {
    stop('`at_ages` must contain the same time slices as `ages`.')
  }
  # find ages if not provided
  if (is.null(ages)) ages <- as.numeric(names(age_list))
  # force ages to be in decreasing order (important in downstream functions)
  ages <- sort(ages, decreasing = TRUE)

  net_list <- list()

  for (m in seq_along(ages)) {
    mat <- age_list[[m]]  # mat is an array when state = 1:2
    mat[mat < threshold] <- 0
    if (!weighted) {
      mat[mat >= threshold] <- 1
    }

    if (length(state) > 1) {
      mat <- apply(mat * slice.index(mat, 3), 1:2, sum)
    }
    mat <- mat[rowSums(mat) != 0, ]
    mat <- mat[, colSums(mat) != 0]
    net_list[[m]] <- mat
  }

  names(net_list) <- ages

  return(net_list)
}
