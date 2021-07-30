#' Find modules in networks sampled across MCMC at specific time slices
#'
#' @param samples_at_ages List of sampled networks at time slices produces by calling posterior_at_ages().
#' @param ages Vector of ages (time slices in the past) at which samples were retrieved.
#'
#' @return Data frame with module membership information for each sampled network at each time slice.
#' @importFrom dplyr mutate group_by distinct summarize left_join case_when select bind_rows tibble
#' @importFrom bipartite empty
#' @export
#'
#' @examples
#' \dontrun{
#'   mod_samples <- modules_from_samples(samples_at_ages, ages)
#' }
modules_from_samples <- function(samples_at_ages, ages) {

  if (length(samples_at_ages) != length(ages)) {
    stop('`samples_at_ages` must contain the same time slices as `ages`.')
  }

  Qsamples <- tibble()
  mod_samples <- tibble()

  nsamp <- dim(samples_at_ages[[1]])[1]

  for(a in 1:length(ages)){
    for(i in 1:nsamp){
      net <- samples_at_ages[[a]][i,,]

      if(ncol(empty(net)) > 1){

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

          mod_samples <- bind_rows(mod_samples, mtbl)
        }
      }
    }
  }

  # remove replicate modules in fully connected networks with identical modules
  duplicates_removed <- remove_duplicate_modules(mod_samples)

  duplicates_removed

}


remove_duplicate_modules <- function(mod_samples) {

  duplicates_removed <- mod_samples %>% group_by(age, sample) %>% distinct(name) %>% summarize(u=n()) %>%
    left_join(mod_samples %>% group_by(age, sample) %>% summarize(n=n())) %>%
    mutate(problem = case_when(u != n ~ "YES", u == n ~ "NO")) %>%
    left_join(mod_samples) %>%
    mutate(original_module = case_when(problem == "YES" ~ 1,
                                       problem == "NO" ~ as.numeric(original_module))) %>%
    distinct() %>%
    select(name, age, sample, original_module)

  duplicates_removed

}

n <- name <- NULL

