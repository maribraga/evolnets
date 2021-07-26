
#' Calculate support for modules from summary networks based on modules of sampled networks
#'
#' Functions to validate modules from summary networks by calculating the frequency with
#' which pairs of nodes are placed in the same module in networks sampled across MCMC
#'
#' @param ages Vector of ages (time points in the past) at which samples will be
#'   retrieved.
#' @param mod_samples Data frame produced by modules_from_samples() containing module membership of each node for each sampled network at each time slice.
#' @param threshold Minimum frequency which two nodes are placed in the same module to consider it supported.
#' @param edge_list Logical. Whether to return a list of edge lists or a list of matrices.
#' @param include_all Logical. Include all nodes or only those present at the time slice?
#' @param modules_across_ages Data frame containing the module information for the summary network.
#'
#' @return list of data frames
#' @importFrom dplyr mutate case_when arrange pull arrange as_tibble bind_rows distinct filter left_join select summarize tibble
#' @importFrom tidyr complete
#' @export
#'
#' @examples
support_for_modules <- function(mod_samples, ages, modules_across_ages, threshold = 0.9, edge_list = TRUE, include_all = FALSE) {

  # calculate pairwise module membership
  pair_mod_tbl <- pairwise_membership(mod_samples, ages, edge_list)

  # find modules with strong support


  # make heatmap
  pair_heatmaps <- list()

  for(i in 1:length(ages)){

    edge_list <- as_tibble(pair_mod_tbl[[i]]) %>%
      purrr::when(include_all ~
                    full_join(., modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
                    full_join(modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name")),
                  ~ inner_join(., modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
                    inner_join(modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name"))) %>%
      mutate(Module = ifelse(module.x == module.y, module.x, NA),
             supported_mod = case_when(freq >= threshold ~ Module))

    order <- edge_list %>% arrange(module.x) %>% pull(row) %>% unique()

    for_heatmap <- edge_list %>% mutate(
      row = factor(row, levels = order),
      col = factor(col, levels = order))

    pair_heatmaps[[i]] <- for_heatmap
  }

  pair_heatmaps

}


pairwise_membership <- function(mod_samples, ages, edge_list = TRUE) {

  pair_mod_matrix <- list() # for heatmap from a matrix
  pair_mod_tbl <- list()    # for heatmap from an edge list

  for(a in 1:length(ages)){
    taxa = mod_samples %>% filter(age == ages[a]) %>% distinct(name)
    ntaxa = nrow(taxa)

    heat <- matrix(data=0, nrow = ntaxa, ncol = ntaxa)
    rownames(heat) <- colnames(heat) <- taxa$name

    tbl <- tibble(row = taxa$name, col = taxa$name) %>%
      complete(row,col) %>%
      mutate(freq = 0)

    mod_samples_at_age <- mod_samples %>% filter(age == ages[a]) %>% distinct(sample) %>% pull(sample)

    for(i in mod_samples_at_age){

      table <- mod_samples %>% filter(age == ages[a], sample == i)
      mods <- unique(table$original_module)

      for(m in mods){
        module <- filter(table, original_module == m)

        for(r in 1:nrow(heat)){
          for(c in 1:ncol(heat)){
            if(rownames(heat)[r] %in% module$name & colnames(heat)[c] %in% module$name){
              heat[r,c] = heat[r,c] + 1
              tbl <- mutate(tbl, freq = case_when(row == rownames(heat)[r] & col == colnames(heat)[c] ~ freq + 1,
                                                  TRUE ~ freq))
            }
          }
        }
      }
    }

    heat <- heat/length(mod_samples_at_age)
    tbl <- mutate(tbl, freq = freq/length(mod_samples_at_age))

    pair_mod_matrix[[a]] <- heat
    pair_mod_tbl[[a]] <- tbl

  }

  if(edge_list){
    return(pair_mod_tbl)
  } else{
    return(pair_mod_matrix)
  }

}

module.y <- module.x <- freq <- name <- NULL
