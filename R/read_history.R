#' Read history.txt file produced by RevBayes
#'
#' @description Read history.txt file containing the posterior samples of interaction histories
#' @param path_to_hist string path to the .txt file exported from RevBayes
#'
#' @return data frame
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # read history file
#' \dontrun{history_full <- read_history("/path/history.txt")}
read_history <- function(path_to_hist){
  colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
  history = utils::read.table(path_to_hist, sep="\t", header=T, colClasses = colclasses) %>%
    dplyr::mutate(node_index = node_index + 1,
           parent_index = parent_index + 1,
           child1_index = child1_index + 1,
           child2_index = child2_index + 1)
}

node_index <- parent_index <- child1_index <- child2_index <- NULL
