#' Get ancestral networks sampled across MCMC
#'
#' @description Retrieve samples in the history.txt file produced by the analysis of host-repertoire evolution in RevBayes as a host-parasite incidence matrix (binary bipartite network).
#'
#' @param history Data frame with posterior samples of interaction histories.
#' @param ages Vector of ages (time points in the past) at which samples will be retrieved.
#' @param tree Parasite tree
#' @param host_tree Host tree
#' @param state Which state? Default is 2. For analyses using the 3-state model, can take the values 1 (potential host) or 2 (actual host).
#' @param drop_empty Logical. Remove taxa without any interactions?
#'
#' @return A list of arrays of samples x parasites x hosts. The number of samples is the number of iterations in `history`. At each age, all hosts and all extant parasite lineages are included.
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
#' # get samples at age
#' ages <- c(80,70)
#' samples_at_ages <- samples_at_ages(history, ages, tree, host_tree)
samples_at_ages <- function(history, ages, tree, host_tree, state=c(2), drop_empty=T) {

  samp_ages = list()
  for (i in 1:length(ages)) {
    age = ages[i]
    samp_ages[[i]] = make_matrix_samples_at_age( history, age, tree, host_tree, state, drop_empty )
  }

  samp_ages
}


make_matrix_samples_at_age = function(dat, age, tree, host_tree, state, drop_empty) {

    iterations = sort(unique(dat$iteration))
    n_iter = length(iterations)

    # get dimensions
    n_host_tip = length( stringr::str_split( dat$start_state[1], "" )[[1]] )
    n_parasite_lineage = length(unique(dat$node_index))
    n_parasite_tip = (n_parasite_lineage + 1) / 2

    m_names = list( 1:n_iter,
                    c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage)),
                    host_tree$tip.label )

    m = array(0, dim=c(n_iter, n_parasite_lineage, n_host_tip), dimnames=m_names)

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
            m[it_idx, n_idx, s_idx ] = 1
        }
    }

    return(m)
}
