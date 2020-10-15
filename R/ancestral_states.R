#' Functions to calculate the posterior probability of ancestral host
#' repertoires
#'
#' @description Group of functions to calculate the posterior probabilities of
#'   ancestral host repertoires at internal nodes of the parasite tree or at
#'   specific time points in the past.
#'
#' @param history Data frame with posterior samples of interaction histories.
#' @param ages Vector of ages (time points in the past) at which ancestral host
#'   repertoires will be reconstructed.
#' @param state Which state? Default is 2. For analyses using the 3-state model,
#'   can take the values 1 (potential host) or 2 (actual host).
#' @param tree Parasite tree
#' @param host_tree Host tree
#' @param drop_empty Logical. Remove taxa without any interactions?
#' @param nodes Vector of internal nodes for which to calculate the posterior
#'   probability of `state`.
#'
#' @examples
#' # read parasite and host tree
#' data(tree)
#' data(host_tree)
#'
#' # read histories sampled during MCMC
#' data(history)
#'
#' # calculate the posterior probability of all possible interactions
#' # at chosen ages (time points in the past)
#' ages <- c(80,70)
#' pp_at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#'
#' # calculate the posterior probability of host repertoires
#' # at chosen internal nodes of the parasite tree
#' nodes <- c(129:131)
#' pp_at_nodes <- posterior_at_nodes(history, nodes, host_tree)
#'
#' @name ancestral_states
NULL

#' @describeIn ancestral_states Make a matrix with marginal posterior probabilities of interactions at given internal nodes.
#' @export
posterior_at_nodes = function(history, nodes, host_tree, state = c(2)) {

  dat <- dplyr::filter(history, node_index %in% nodes)
  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)

  # get dimensions
  n_host_tip = length( stringr::str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  g = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)

  for (it in iterations) {
    dat_it = dat[ dat$iteration == it, ]
    ret = list()
    for (i in 1:length(nodes)) {
      dat2 = dat_it[ dat_it$node_index == nodes[i], ]
      if(nrow(dat2)==1){
        ret[[i]] = dat2
      } else{
        ret[[i]] = dat2[ which.min(dat2$transition_time), ]
      }
    }

    ret <- data.table::rbindlist(ret)

    for (r in 1:nrow(ret)) {
      s = as.numeric( stringr::str_split (ret$end_state[r], "")[[1]] )
      s_idx = s %in% state

      g[ r, s_idx ] = g[ r, s_idx ] + 1
    }
  }

  # convert to probability
  g = g * (1/n_iter)
  row.names(g) <- paste0("Index_",nodes)
  colnames(g) <- host_tree$tip.label

  return(g)
}


#' @describeIn ancestral_states Make a list of matrices with posterior probabilities for each interaction of extant lineages at given ages.
#' @export
posterior_at_ages <- function(history, ages, tree, host_tree, state=c(2), drop_empty=T) {

  list_m_at_ages = list()

  for (i in 1:length(ages)) {
    age = ages[i]
    list_m_at_ages[[i]] = t(make_matrix_at_age( history, age, tree, host_tree, state, drop_empty ))
  }

  list_m_at_ages
}


# Make matrix with marginal posterior probabilities of interactions at a given age
make_matrix_at_age = function(dat, age, tree, host_tree, state, drop_empty) {

  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)

  # get dimensions
  n_host_tip = length( stringr::str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  n_parasite_tip = (n_parasite_lineage + 1) / 2
  m = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)

    for (it in iterations) {

        # get dataset for iteration
        dat_it = dat[ dat$iteration==it, ]

        # extract relevant branches
        dat_it = make_dat_age(dat_it, age)

        # add edges ( parasite x host )
        for (i in 1:nrow(dat_it)) {
            n_idx = dat_it$node_index[i]
            s = as.numeric( stringr::str_split (dat_it$end_state[i], "")[[1]] )
            s_idx = s %in% state
            m[ n_idx, s_idx ] = m[ n_idx, s_idx ] + 1
        }
    }

    row.names(m) <- c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage))
    colnames(m) <- host_tree$tip.label

    # remove empty rows/columns
    if (drop_empty) {
        m = m[ rowSums(m)!=0, ]
        m = m[ ,colSums(m)!=0 ]
    }

    # convert to probability
    m = m * (1/n_iter)

    return(m)
}


# find lineages (and their repertoires) that exist during the specified age
make_dat_age = function(dat, age) {

    # reduce the history data frame to only include relevant branches
    dat2 = dat[dat$branch_start_time>=age & dat$branch_end_time<=age, ]

    nodes = unique(dat2$node_index)

    dat3 = rbind( dat2[ dat2$transition_type=="no_change", ],
                  dat2[ dat2$transition_type=="anagenetic" & dat2$transition_time>=age, ] )

    ret = list()
    for (i in 1:length(nodes)) {
      if(!(nodes[i] %in% dat3$node_index)){               # if changes happened only after age
        parent = dat[ dat$node_index==nodes[i], 12][1]    # get state at parent node
        dat4 = dat[ dat$node_index==parent, ]
        dat4$node_index = nodes[i]                        # fix node_index - back to child nodes
        if(nrow(dat4)==1){                                # when type is no_change, time is NA and which.min doesn't work
          ret[[i]] = dat4
        } else{
          ret[[i]] = dat4[ which.min(dat4$transition_time), ]
        }
      } else {
      dat4 = dat3[ dat3$node_index==nodes[i], ]
      if(nrow(dat4)==1){
        ret[[i]] = dat4
      } else{
      ret[[i]] = dat4[ which.min(dat4$transition_time), ]  # get state at minimum, not maximum time (that is greater than age)
      }
      }
    }
    return(data.table::rbindlist(ret))
}


