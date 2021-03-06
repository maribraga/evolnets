#' Calculate the rate of host-repertoire evolution
#'
#' Calculate the number of host gains, host losses, and the
#'   effective rate of host repertoire evolution
#' @param history A data frame containing the character history produced by
#'   RevBayes
#' @param tree A phylogenetic tree of the parasite clade
#'
#' @examples
#' # read parasite tree
#' data(tree)
#'
#' # read hist
#' data(history)
#'
#' # remove burnin and thin samples
#' it_seq <- seq(10000,100000,1000)
#' history <- dplyr::filter(history, iteration %in% it_seq)
#'
#' # effective rate
#' n_events <- count_events(history)
#' rate <- effective_rate(history,tree)
#'
#' # gains and losses
#' gl_events <- count_gl(history)
#' gl_rates <- rate_gl(history, tree)
#'
#' @name events_counter
NULL

#' @describeIn events_counter Get the total number of events along the parasite tree
#' @export
count_events <- function(history){

  root <- max(history$node_index)

  n_events <- history %>%
    dplyr::filter(transition_type == "anagenetic", node_index < root) %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(n = dplyr::n(), .groups = 'drop') %>%
    dplyr::summarise(mean = mean(n)) %>%
    dplyr::pull(mean)

  n_events

}


#' @describeIn events_counter Get the effective rate of evolution, i.e. number
#'   of events per million years, along each tree branch
#' @export
effective_rate <- function(history, tree) {

  nev <- count_events(history)
  tree_length <- sum(tree$edge.length)
  rate <- nev/tree_length

  rate
}


#' @describeIn events_counter Get the number of host gains and host losses
#' @export
count_gl <- function(history) {

  root <- max(history$node_index)

  ngl <- history %>%
    dplyr::filter(transition_type == "anagenetic", node_index < root) %>%
    dplyr::select(iteration, start_state, end_state) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(start = sum(as.numeric(stringr::str_split(start_state, "")[[1]])),
                  end = sum(as.numeric(stringr::str_split(end_state, "")[[1]]))) %>%
    dplyr::mutate(type = dplyr::case_when(end > start ~ "gain",
                                          end < start ~ "loss"))

  ngains <- dplyr::filter(ngl, type == "gain") %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::summarise(mean = mean(n)) %>%
    dplyr::pull(mean)

  nloss <- dplyr::filter(ngl, type == "loss") %>%
    dplyr::group_by(iteration) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::summarise(mean = mean(n)) %>%
    dplyr::pull(mean)

  gl <- c(ngains, nloss)
  names(gl) <- c("gains","losses")

  gl

}


#' @describeIn events_counter Get the effective rate of host gain and host loss
#' @export
rate_gl <- function(history, tree) {

  gl <- count_gl(history)
  tree_length <- sum(tree$edge.length)
  rg <- gl[1]/tree_length
  rl <- gl[2]/tree_length

  rgl <- c(rg,rl)
  names(rgl) <- c("gain_rate","loss_rate")

  rgl

}


iteration <- transition_type <- node_index <- n <- start_state <- end_state <- type <- NULL
