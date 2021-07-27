#' Functions to calculate the posterior probability of ancestral host
#' repertoires at internal nodes
#'
#' Group of functions to calculate the posterior probabilities of
#'   ancestral host repertoires at internal nodes of the parasite tree.
#'
#' @param history Data frame with posterior samples of interaction histories.
#' @param state Which state? Default is 2. For analyses using the 3-state model,
#'   can take the values 1 (potential host) or 2 (actual host).
#' @param host_tree Host tree
#' @param nodes Vector of internal nodes for which to calculate the posterior
#'   probability of `state`.
#'
#' @return A matrix with marginal posterior probabilities of interactions at given internal nodes.
#' @export
#'
#' @examples
#' # read parasite and host tree
#' data(tree)
#' data(host_tree)
#'
#' # read histories sampled during MCMC
#' data(history)
#'
#' # calculate the posterior probability of host repertoires
#' # at chosen internal nodes of the parasite tree
#' nodes <- c(129:131)
#' pp_at_nodes <- posterior_at_nodes(history, nodes, host_tree)
posterior_at_nodes <- function(history, nodes, host_tree, state = c(2)) {
  dat <- dplyr::filter(history, node_index %in% nodes)
  iterations <- sort(unique(dat$iteration))
  n_iter <- length(iterations)

  # get dimensions
  n_host_tip <- length(stringr::str_split( dat$start_state[1], "" )[[1]])
  n_parasite_lineage <- length(unique(dat$node_index))

  g <- matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)

  array_names <- list(1:n_iter, paste0("Index_",nodes), host_tree$tip.label)
  array <- array(0, dim = c(n_iter, n_parasite_lineage, n_host_tip), dimnames = array_names)

  for (i in 1:n_iter) {
    it <- iterations[i]
    dat_it <- dat[dat$iteration == it, ]
    ret <- list()
    for (i in 1:length(nodes)) {
      dat2 <- dat_it[dat_it$node_index == nodes[i], ]
      if (nrow(dat2) == 1) {
        ret[[i]] <- dat2
      } else {
        ret[[i]] <- dat2[which.min(dat2$transition_time), ]
      }
    }

    ret <- data.table::rbindlist(ret)

    for (r in 1:nrow(ret)) {
      s <- as.numeric(stringr::str_split(ret$end_state[r], "")[[1]])
      s_idx = s %in% state

      g[r, s_idx ] <- g[r, s_idx ] + 1
      array[i, r, s_idx] <- state
    }
  }

  # convert to probability
  g = g * (1 / n_iter)
  row.names(g) <- paste0("Index_", nodes)
  colnames(g) <- host_tree$tip.label

  list <- list(array, g)

  return(list)
}

