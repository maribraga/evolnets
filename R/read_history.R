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
    dplyr::filter(history, iteration > max(iteration)*burnin)
}

node_index <- parent_index <- child1_index <- child2_index <- NULL
