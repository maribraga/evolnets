
#' Validate modules from summary networks
#'
#' Functions to validate modules from summary networks by calculating the frequency with
#' which pairs of nodes are placed in the same module in networks sampled across MCMC
#'
#' @param samples_at_ages
#' @param ages
#' @param modules_across_ages
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
validate_modules <- function(samples_at_ages, ages, modules_across_ages, include_all = FALSE, seed = NULL, edge_list = TRUE) {

  # find modules for each sample at each age
  Mod_samples <- modules_from_samples(samples_at_ages, ages, seed)

  # remove fully connected networks with identical modules
  good_samples <- remove_bad_module_samples(Mod_samples)

  # calculate pairwise module membership
  pair_mod_tbl <- pairwise_membership(Mod_samples, ages, good_samples, edge_list)

  # match modules across ages
  pair_heatmaps <- list()

  for(i in 1:length(ages)){

    edge_list <- as_tibble(pair_mod_tbl[[i]]) %>%
      purrr::when(include_all ~
                    full_join(., modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
                    full_join(modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name")),
                  ~ inner_join(., modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
                    inner_join(modules_across_ages %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name"))) %>%
      mutate(Module = ifelse(module.x == module.y, module.x, NA), Frequency = freq/max(freq))

    order <- edge_list %>% arrange(module.x) %>% pull(row) %>% unique()

    for_heatmap <- edge_list %>% mutate(
      row = factor(row, levels = order),
      col = factor(col, levels = order))

    pair_heatmaps[[i]] <- for_heatmap
  }

  pair_heatmaps

}

modules_from_samples <- function(samples_at_ages, ages, seed = NULL) {

  Qsamples <- tibble()
  Mod_samples <- tibble()

  nsamp <- dim(samples_at_ages[[1]])[1]

  for(a in 1:length(ages)){
    for(i in 1:nsamp){
      net <- samples_at_ages[[a]][i,,]

      if(ncol(empty(net)) > 1){

        set.seed(seed)
        mod <- mycomputeModules(net)

        q <- mod@likelihood
        Qsamples <- bind_rows(Qsamples, tibble(age = ages[a], sample = i, Q=q))

        mod_list <- listModuleInformation(mod)[[2]]
        nmod <- length(mod_list)
        for(m in 1:nmod){
          members <- unlist(mod_list[[m]])
          mtbl <- tibble(name = members,
                         age = rep(ages[a], length(members)),
                         sample = rep(i, length(members)),
                         original_module = rep(m, length(members)))

          Mod_samples <- bind_rows(Mod_samples, mtbl)
        }
      }
    }
  }

  Mod_samples

}

#modules_across_ages <- function()

remove_bad_module_samples <- function(Mod_samples) {

  good_samples <- Mod_samples %>% group_by(age, sample) %>% distinct(name) %>% summarize(u=n()) %>%
    left_join(Mod_samples %>% group_by(age, sample) %>% summarize(n=n())) %>%
    mutate(problem = case_when(u != n ~ "YES", u == n ~ "NO")) %>%
    filter(problem == "NO") %>%
    select(age,sample)

  good_samples

}

pairwise_membership <- function(Mod_samples, ages, good_samples, edge_list = TRUE) {

  pair_mod_matrix <- list() # for heatmap from a matrix
  pair_mod_tbl <- list()    # for heatmap from an edge list

  for(a in 1:length(ages)){
    taxa = Mod_samples %>% filter(age == ages[a]) %>% distinct(name)
    ntaxa = nrow(taxa)

    heat <- matrix(data=0, nrow = ntaxa, ncol = ntaxa)
    rownames(heat) <- colnames(heat) <- taxa$name

    tbl <- tibble(row = taxa$name, col = taxa$name) %>%
      complete(row,col) %>%
      mutate(freq = 0)

    good_samples_at_age <- good_samples %>% filter(age == ages[a]) %>% pull(sample)

    for(i in good_samples_at_age){

      table <- Mod_samples %>% filter(age == ages[a], sample == i)
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

  pair_mod_matrix[[a]] <- heat
  pair_mod_tbl[[a]] <- tbl

  }

  if(edge_list){
    return(pair_mod_tbl)
  } else{
    return(pair_mod_matrix)
  }

}

