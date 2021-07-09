
#' Title
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
validate_modules <- function(samples_at_ages, ages, modules_across_ages, seed = NULL) {

  Mod_samples <- modules_from_samples(samples_at_ages, ages, modules_across_ages, seed)

  good_samples <- remove_bad_module_samples(Mod_samples)

  pair_mod_tbl <- pairwise_membership(Mod_samples, ages, good_samples)

  pair_heatmap

}

modules_from_samples <- function(samples_at_ages, ages, modules_across_ages, seed = NULL) {

  Qsamples <- tibble()
  Mod_samples <- tibble()


  Mod_samples

}


remove_bad_module_samples <- function(Mod_samples) {

  out <- Mod_samples %>% group_by(age, sample) %>% distinct(name) %>% summarize(u=n()) %>%
  left_join(Mod_samples %>% group_by(age, sample) %>% summarize(n=n())) %>%
  mutate(problem = case_when(u != n ~ "YES", u == n ~ "NO")) %>%
  filter(problem == "YES")


  good_samples

}

pairwise_membership <- function(Mod_samples, ages, good_samples) {

  pair_mod_matrix <- list() # for heatmap from a matrix
  pair_mod_tbl <- list()    # for heatmap from an edge list

  pair_mod_tbl

}

