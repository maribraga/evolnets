#' Read history.txt file produced by RevBayes
#'
#' Read history.txt file containing the posterior samples of
#'   interaction histories
#'
#' @param path_to_hist String path to the .txt file exported from RevBayes
#' @param burnin Fraction of iterations to be removed as burnin. Default to 10%.
#'
#' @return A data.frame
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # read history file
#' \dontrun{history <- read_history("/path/history.txt", burnin = 0.2)}

read_history <- function(path_to_hist, burnin = 0.1){

  if (!is.numeric(burnin) || burnin > 1 || burnin < 0 || length(burnin) != 1) {
    stop('`burnin` should be a numeric vector of length 1, between 0 and 1.')
  }

  colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))

  history <- utils::read.table(path_to_hist, sep = "\t", header = TRUE, colClasses = colclasses) %>%
    dplyr::mutate(
      node_index = .data$node_index + 1,
      parent_index = .data$parent_index + 1,
      child1_index = .data$child1_index + 1,
      child2_index = .data$child2_index + 1
    ) %>%
    dplyr::filter(.data$iteration > max(.data$iteration) * burnin)
}


#' Read symbiont tree file produced by RevBayes containing node indices
#'
#' @param path_to_tree String path to the .txt file exported from RevBayes
#'
#' @return A `phylo` object with node names given by RevBayes, which are important to place the inferred ancestral states in the symbiont tree.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{tree <- read_tree_from_revbayes("symbiont_tree.tre")}

read_tree_from_revbayes <- function(path_to_tree){

  treeRev <- treeio::read.beast.newick(path_to_tree)

  tree <- treeRev@phylo

  index_node <- treeRev@data %>%
    dplyr::mutate(node = as.numeric(.data$node)) %>%
    dplyr::arrange(.data$node)

  indices <- index_node %>%
    dplyr::filter(.data$node > treeio::Ntip(tree)) %>%
    dplyr::pull(.data$index)

  names(indices) <- NULL

  tree$node.label <- paste0("Index_", indices)

  return(tree)
}
