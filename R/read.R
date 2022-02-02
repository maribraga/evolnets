#' Read history.txt file produced by RevBayes
#'
#' Read history.txt file containing the posterior samples of
#'   interaction histories
#'
#' @param path_to_hist String path to the .txt file exported from RevBayes
#' @param burnin Fraction of iterations to be removed as burnin. Default to 10%.
#'
#' @return A data frame
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # read history file
#' \dontrun{history <- read_history("/path/history.txt", burnin = 0.2)}

read_history <- function(path_to_hist, burnin = 0.1){
  colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))

  history = utils::read.table(path_to_hist, sep="\t", header=T, colClasses = colclasses) %>%
    dplyr::mutate(node_index = node_index + 1,
                  parent_index = parent_index + 1,
                  child1_index = child1_index + 1,
                  child2_index = child2_index + 1) %>%
    dplyr::filter(iteration > max(iteration)*burnin)
}

node_index <- parent_index <- child1_index <- child2_index <- iteration <- NULL


#' Read symbiont tree file produced by RevBayes containing node indices
#'
#' @param path_to_tree String path to the .txt file exported from RevBayes
#'
#' @return A `phylo` object
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{tree <- read_tree_from_revbayes("symbiont_tree.tre")}

read_tree_from_revbayes <- function(path_to_tree){

  treeRev <- treeio::read.beast.newick(path_to_tree)

  tree <- treeRev@phylo

  index_node <- treeRev@data %>%
    dplyr::mutate(node = as.numeric(node)) %>%
    dplyr::arrange(node)

  indices <- index_node %>%
    dplyr::filter(node > treeio::Ntip(tree)) %>%
    dplyr::pull(index)

  names(indices) <- NULL

  tree$node.label <- paste0("Index_",indices)

  return(tree)
}
