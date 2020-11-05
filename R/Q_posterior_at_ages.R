#' Posterior distribution of modularity across ancestral networks
#'
#' Calculate z-scores for modularity (Q) for each MCMC sample at
#'   time points in the past based on null networks where all interactions have
#'   the same probability. By calculating z-scores, we can compare ancestral
#'   networks at different ages.
#'
#' @param samples_at_ages List of ancestral networks sampled across MCMC at
#'   given ages.
#' @param ages Vector of ages (time points in the past) of ancestral networks.
#' @param null Number of null networks to generate to calculate the z-score.
#' @param seed Seed passed to `stats::simulate` to generate null networks.
#'   Default to NULL.
#'
#' @return A tibble of Q z-scores and p-values across samples and ages.
#' @importFrom magrittr %>%
#' @importFrom methods is slot
#' @export
#'
#' @examples
Q_posterior_at_ages <- function(samples_at_ages, ages, null = 100, seed = NULL){

  Q_null <- Q_samples_null(samples_at_ages, ages, null, seed)
  Q_samples <- Q_samples_at_ages(samples_at_ages, ages)

  Qzsamples <- Q_null %>%
    dplyr::group_by(age, sample) %>%
    dplyr::summarize(mean_Q = mean(Qnull), sd_Q = stats::sd(Qnull), .groups = 'drop') %>%
    dplyr::left_join(Q_samples, by = c("age", "sample")) %>%
    dplyr::mutate(z = (obs_Q - mean_Q)/sd_Q) %>%
    dplyr::left_join(Q_null %>%
                       dplyr::left_join(Q_samples, by = c("age", "sample")) %>%
                       dplyr::group_by(age, sample) %>%
                       dplyr::summarise(p = sum(Qnull >= obs_Q)/null, .groups = 'drop'),
                     by = c("age", "sample"))

  Qzsamples

}


# Simulate null networks and calculate Q
Q_samples_null <- function(samples_at_ages, ages, null = 100, seed = NULL){

  nnull <- null
  nsamp <- dim(samples_at_ages[[1]])[1]
  Q_null <- dplyr::tibble()

  for(a in 1:length(samples_at_ages)){

    Nulls_age <- list()

    for(i in 1:nsamp){
      net <- samples_at_ages[[a]][i,,]
      net = net[ rowSums(net)!=0, ]
      net = net[ ,colSums(net)!=0 ]

      nullm <- vegan::nullmodel(net, "r00")
      sim <- stats::simulate(nullm, nsim=nnull, seed = seed)
      Nulls_age[[i]] <- sim

      for(j in 1:nnull){
        mod <- mycomputeModules(sim[,,j])
        Qrandom <- mod@likelihood

        Q_null <- dplyr::bind_rows(Q_null, dplyr::tibble(age = ages[a],
                                                               sample = i,
                                                               sim = j,
                                                               Qnull = Qrandom))
      }
    }
  }

  Q_null

}


# Get Q for each sampled network
Q_samples_at_ages <- function(samples_at_ages, ages){

  nsamp <- dim(samples_at_ages[[1]])[1]
  Q_samples <- dplyr::tibble()

  for(a in 1:length(samples_at_ages)){
    for(i in 1:nsamp){
      net <- samples_at_ages[[a]][i,,]
      mod <- mycomputeModules(net)
      Q <- mod@likelihood

      Q_samples <- dplyr::bind_rows(Q_samples, dplyr::tibble(age = ages[a], sample = i, obs_Q = Q))
    }
  }

  Q_samples

}


age <- sample <- obs_Q <- Qnull <- mean_Q <- sd_Q <- NULL
