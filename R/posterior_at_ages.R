#' Functions to calculate the posterior probability of ancestral host
#' repertoires at ages.
#'
#' Group of functions to retrieve samples in the history.txt file produced by the
#'   analysis of host-repertoire evolution in RevBayes and calculate the posterior
#'   probabilities of ancestral host repertoires at specific time points in the past.
#'
#'
#' @param history Data frame with posterior samples of interaction histories.
#' @param ages Vector of ages (time points in the past) at which samples will be
#'   retrieved.
#' @param tree Parasite tree
#' @param host_tree Host tree
#' @param state Which state? Default is 2. For analyses using the 3-state model,
#'   can take the values 1 (potential host) or 2 (actual host).
#' @param drop_empty Logical. Remove taxa without any interactions?
#'
#' @return A list of arrays of samples x parasites x hosts (first element) and a matrix of
#'   posterior probabilities (second element). The number of samples is the number of
#'   iterations in `history`. At each age, all hosts and all extant parasite lineages are included.
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
#' # get samples at ages
#' ages <- c(80,70)
#' at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#' samples_at_ages <- at_ages[[1]]
#'
#' # get posterior probabilities at ages
#' pp_at_ages <- at_ages[[2]]
posterior_at_ages <- function(history, ages, tree, host_tree, state=c(2), drop_empty=T) {

  ages_samp_post = list()
  samp_ages = list()
  post_ages = list()
  samp_post_ages = list()

  for (i in 1:length(ages)) {
    age = ages[i]
    ages_samp_post = make_samples_post_at_age( history, age, tree, host_tree, state, drop_empty )

    samp_ages[[i]] = ages_samp_post[[1]]
    post_ages[[i]] = ages_samp_post[[2]]
  }

  samp_post_ages <- list(samp_ages, post_ages)
}


make_samples_post_at_age = function(dat, age, tree, host_tree, state, drop_empty) {

  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)

  # get dimensions
  n_host_tip = length( stringr::str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  n_parasite_tip = (n_parasite_lineage + 1) / 2

  m_names = list( 1:n_iter,
                  c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage)),
                  host_tree$tip.label )

  m_sample = array(0, dim=c(n_iter, n_parasite_lineage, n_host_tip), dimnames=m_names)

  m_posterior = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)

  for (it_idx in 1:length(iterations)) {
    it = iterations[it_idx]

    # get dataset for iteration
    dat_it = dat[ dat$iteration==it, ]

    # extract relevant branches
    dat_it = make_dat_age(dat_it, age)

    # add edges ( parasite x host )
    for (i in 1:nrow(dat_it)) {
      n_idx = dat_it$node_index[i]
      s = as.numeric( stringr::str_split (dat_it$end_state[i], "")[[1]] )
      s_idx = s %in% state

      m_sample[it_idx, n_idx, s_idx ] = state

      m_posterior[ n_idx, s_idx ] = m_posterior[ n_idx, s_idx ] + 1
    }
  }

  row.names(m_posterior) <- c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage))
  colnames(m_posterior) <- host_tree$tip.label

  # remove empty rows/columns
  if (drop_empty) {
    m_posterior = m_posterior[ rowSums(m_posterior)!=0, ]
    m_posterior = m_posterior[ ,colSums(m_posterior)!=0 ]
  }

  # convert to probability
  m_posterior = m_posterior * (1/n_iter)

  samp_post = list(m_sample, m_posterior)

  return(samp_post)
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

