#' Read history.txt file produced by RevBayes
#'
#' @param path_to_hist string
#'
#' @return data frame
#' @export
#'
#' @examples
read_history <- function(path_to_hist){
  colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
  history = utils::read.table(path_to_hist, sep="\t", header=T, colClasses = colclasses) %>%
    dplyr::mutate(node_index = node_index + 1,
           parent_index = parent_index + 1,
           child1_index = child1_index + 1,
           child2_index = child2_index + 1)
}
