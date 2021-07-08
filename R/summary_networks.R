#' Get summary networks from posterior probabilities of interactions of extant
#' lineages at given time points in the past
#'
#' @param pp_at_ages List of matrices with posterior probabilities for each
#'   interaction of extant lineages at given ages
#' @param pt Probability threshold to include an interaction in the network.
#'   Interactions with pp < pt will be dropped.
#' @param weighted Logical. Use posterior probabilities as interaction weights?
#'
#' @return summary network (incidence matrix)
#' @export
#'
#' @examples
#' #' data(tree)
#' data(host_tree)
#' data(history)
#'
#' ages <- c(60,50,40)
#' pp_at_ages <- posterior_at_ages(history, ages, tree, host_tree)[[2]]
#' weighted_net_50 <- get_summary_network(pp_at_ages, pt = 0.5, weighted = TRUE)
#' binary_net_90 <- get_summary_network(pp_at_ages, pt = 0.9, weighted = FALSE)
get_summary_network <- function(pp_at_ages, pt, weighted = TRUE){

  net_list <- list()

  for(m in 1:length(pp_at_ages)){
    matrix <- pp_at_ages[[m]]
    for(i in 1:nrow(matrix)){
      for(j in 1:ncol(matrix)){
        if(matrix[i,j] < pt){
          matrix[i,j] = 0
        } else{
          if(weighted == FALSE){
            matrix[i,j] = 1
          }
        }
      }
    }
    df <- as.data.frame(matrix)
    df = df[ rowSums(df)!=0, ]
    df = df[ ,colSums(df)!=0 ]
    net_list[[m]] <- df
  }

  net_list
}


