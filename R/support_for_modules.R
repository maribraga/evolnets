
#' Calculate support for modules from summary networks based on modules of sampled networks
#'
#' Validate modules from summary networks by calculating the frequency with
#' which pairs of nodes are placed in the same module in networks sampled across MCMC
#'
#' @param mod_samples Data frame produced by `modules_from_samples()` containing module membership of each node for each sampled network at each time slice.
#' @param ages Vector of ages (time points in the past) at which samples will be retrieved.
#' @param threshold Minimum frequency with which two nodes are placed in the same module to consider it supported.
#' @param edge_list Logical. Whether to return a list of edge lists or a list of matrices of pairwise frequency.
#' @param include_all Logical. Include all nodes or only those present at the time slice?
#' @param modules_across_ages Data frame containing the module information for the summary network.
#' @param palette Optional. Color palette for module colors in the plot.
#'
#' @return A list containing:
#' 1) `plot`: A plot of pairwise frequency with which the two nodes are placed in the same module at each time slice in `ages`;
#' 2) `pairwise_membership`: A list of edge lists or matrices with the pairwise frequencies at each time slice;
#' 3) `mean_support`: A list of mean and geometric mean pairwise frequency for each module at each time slice.
#' @importFrom dplyr mutate case_when arrange pull arrange as_tibble bind_rows distinct filter left_join select summarize tibble desc
#' @importFrom tidyr complete
#' @importFrom rlang .data
#' @importFrom patchwork wrap_plots
#' @importFrom stats reorder
#' @export
#'
#' @examples
#' \dontrun{
#'   # find modules for each sampled network
#'   mod_samples <- modules_from_samples(samples_at_ages, ages)
#'
#'   # calculate support
#'   support <- support_for_modules(mod_samples, ages, modules_across_ages)
#'   support$plot
#' }
support_for_modules <- function(mod_samples, ages, modules_across_ages, threshold = 0.7, edge_list = TRUE, include_all = FALSE, palette = NULL) {

  if (length(unique(mod_samples$age)) != length(ages)) {
    stop('`mod_samples` must contain the same time slices as `ages`.')
  }
  if (length(unique(modules_across_ages$age)) < length(ages)) {
    stop('`modules_across_ages` must contain all time slices in `ages`.')
  }

  # calculate pairwise module membership
  pair_mod_tbl <- pairwise_membership(mod_samples, ages, edge_list)

  # make heatmaps
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

  names(pair_heatmaps) <- ages

  # calculate mean support
  means <- list()

  for(i in 1:length(ages)){
    means_age <- data.frame(module = NULL, mean = NULL, geo_mean = NULL)
    mods <- sort(unique(pair_heatmaps[[i]]$module.x))

    for(m in mods){
      within_mod <- filter(pair_heatmaps[[i]], .data$Module == m)
      mean <- mean(within_mod$freq)
      gmean <- exp(mean(log(within_mod$freq)))

      means_age <- bind_rows(means_age, data.frame(module = m, mean = mean, geo_mean = gmean))
    }

    means[[i]] <- means_age
  }

  names(means) <- ages

  plot <- plot_pairwise_membership(pair_heatmaps, ages)

  support_list <- list(plot, pair_heatmaps, means)
  names(support_list) <- c("plot", "pairwise_membership", "mean_support")

  support_list

}


plot_pairwise_membership <- function(pair_heatmaps, ages, palette = NULL){

  nages <- length(pair_heatmaps)

  for(a in 1:nages){
    heatmap <- pair_heatmaps[[a]]

    p <- ggplot2::ggplot(heatmap, ggplot2::aes(x = row, y = reorder(col,desc(col)),
                             fill = .data$supported_mod,
                             alpha = freq)) +
      ggplot2::geom_tile() +
      ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::scale_alpha(ggplot2::aes(range = c(min(freq), max(freq)))) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 270),
        legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(paste0(ages[a]," Ma"))

    if(!is.null(palette)){
      p <- p + ggplot2::scale_fill_manual(values = palette, na.value = "grey20", drop = F)
    }

    if(a == 1){
      plot <- p
    } else{
      plot <- wrap_plots(plot, p)
    }

  }

  plot

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
        module <- filter(table, .data$original_module == m)

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
