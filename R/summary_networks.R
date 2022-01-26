#' Get summary networks from posterior probabilities of interactions of extant
#' lineages at given time points in the past
#'
#' @param pp_at_ages List of matrices with posterior probabilities for each interaction of extant
#'   lineages at given ages. Usually calculated `posterior_at_ages`, see example.
#' @param pt Probability threshold to include an interaction in the network.
#'   Interactions with posterior `probability < pt` will be dropped.
#' @param ages Vector of ages (time points in the past) at which samples were retrieved. By default,
#'   uses all ages present in `pp_at_ages`.
#' @param weighted Logical. Use posterior probabilities as interaction weights?
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
#' pp_at_ages <- at_ages$posterior_probabilities
#' weighted_net_50 <- get_summary_network(pp_at_ages, pt = 0.5, weighted = TRUE)
#' binary_net_90 <- get_summary_network(pp_at_ages, pt = 0.9, weighted = FALSE)
get_summary_network <- function(pp_at_ages, pt, ages = NULL, weighted = TRUE){

  if (!is.list(pp_at_ages)) stop('`pp_at_ages` should be a list.')
  if (!is.numeric(pt) || !(pt > 0 & pt <= 1)) stop('`pt` should be a numeric value between 0 and 1.')
  if (!is.null(ages) && (length(pp_at_ages) != length(ages))) {
    stop('`pp_at_ages` must contain the same time slices as `ages`.')
  }
  if (!is.null(ages) && !is.numeric(ages)) stop('`ages` should be a numeric vector or NULL.')
  if (!is.logical(weighted)) stop('`logical` should be a logical value (TRUE/FALSE).')

  # find ages if not provided
  if (is.null(ages)) ages <- as.numeric(names(pp_at_ages))
  net_list <- list()

  for (m in seq_along(ages)) {
    mat <- pp_at_ages[[m]]
    mat[mat < pt] <- 0
    if (!weighted) {
      mat[mat >= pt] <- 1
    }

    df <- as.data.frame(mat)
    df <- df[rowSums(df) != 0, ]
    df <- df[, colSums(df) != 0]
    net_list[[m]] <- df
  }

  names(net_list) <- ages

  net_list
}
