#' Posterior distribution of network structure indices across ancestral networks
#'
#' Calculate z-scores for nestedness (NODF) or modularity (Q) for each MCMC sample at
#'   time points in the past based on null networks where all interactions have
#'   the same probability. By calculating z-scores, we can compare ancestral
#'   networks at different ages.
#'
#' @param samples_at_ages List of ancestral networks sampled across MCMC at
#'   given ages.
#' @param ages Vector of ages (time points in the past) of ancestral networks.
#' @param index Index to be calculated for each ancestral network. "NODF" to
#' calculate nestedness or "Q" to calculate modularity.
#' @param null Number of null networks to generate to calculate the z-score.
#' @param seed Seed passed to `stats::simulate` to generate null networks and set before calculating Q.
#'   Default to NULL.
#'
#' @return A data.frame of z-scores and p-values across samples and ages.
#' @importFrom magrittr %>%
#' @importFrom methods is slot
#' @export
#'
#' @examples
#' data(tree)
#' data(host_tree)
#' data(history)
#'
#' ages <- c(60,50,40)
#' samples_at_ages <- posterior_at_ages(history, ages, tree, host_tree)$Samples
#'
#' # calculate posterior distribution of nestedness
#' Nz <- index_at_ages(samples_at_ages, ages, index = "NODF")
#'
#' #  calculate posterior distribution of modularity (SLOW!)
#' # Qz <- index_at_ages(samples_at_ages, ages, index = "Q")
index_at_ages <- function(samples_at_ages, ages, index, null = 100, seed = NULL){

  if(index == "NODF"){

  NODF_null <- NODF_samples_null(samples_at_ages, ages, null, seed)
  NODF_samples <- NODF_samples_at_ages(samples_at_ages, ages)

  Nzsamples <- NODF_null %>%
    dplyr::group_by(age, sample) %>%
    dplyr::summarize(mean_NODF = mean(NODFnull), sd_NODF = stats::sd(NODFnull), .groups = 'drop') %>%
    dplyr::left_join(NODF_samples, by = c("age", "sample")) %>%
    dplyr::mutate(z = (obs_NODF - mean_NODF)/sd_NODF) %>%
    dplyr::left_join(NODF_null %>%
                       dplyr::left_join(NODF_samples, by = c("age", "sample")) %>%
                       dplyr::group_by(age, sample) %>%
                       dplyr::summarise(p = sum(NODFnull >= obs_NODF)/null, .groups = 'drop'),
                     by = c("age", "sample"))

  Nzsamples

  } else if(index == "Q"){

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

    as.data.frame(Qzsamples)

  } else stop("index must match one of the available indices")

}


# Simulate null networks and calculate NODF
NODF_samples_null <- function(samples_at_ages, ages, null = 100, seed = NULL){

  nnull <- null
  nsamp <- dim(samples_at_ages[[1]])[1]
  NODF_null <- dplyr::tibble()

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
        Nrandom <- bipartite::networklevel(sim[,,j],index="NODF")
        NODF_null <- dplyr::bind_rows(NODF_null, dplyr::tibble(age = ages[a],
                                                               sample = i,
                                                               sim = j,
                                                               NODFnull = Nrandom))
      }
    }
  }

  NODF_null

}


# Get NODF for each sampled network
NODF_samples_at_ages <- function(samples_at_ages, ages){

  nsamp <- dim(samples_at_ages[[1]])[1]
  NODF_samples <- dplyr::tibble()

  for(a in 1:length(samples_at_ages)){
    for(i in 1:nsamp){
      net <- samples_at_ages[[a]][i,,]
      nodf <- bipartite::networklevel(net, index="NODF")
      NODF_samples <- dplyr::bind_rows(NODF_samples, dplyr::tibble(age = ages[a], sample = i, obs_NODF=nodf))
    }
  }

  NODF_samples

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
        set.seed(seed)
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

age <- sample <- obs_NODF <- NODFnull <- mean_NODF <- sd_NODF <- obs_Q <- Qnull <- mean_Q <- sd_Q <- NULL

