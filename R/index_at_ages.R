#' Plot network structure indices across ancestral summary, sampled and extant networks
#'
#' Plot z-scores for nestedness (NODF) and/or modularity (Q) for sampled and summary networks at time points in
#' the past, calculated by `index_at_ages_samples` and `index_at_ages_summary`.
#'
#' @param nodf_sampled Output of `index_at_ages_samples` when index = "NODF".
#' @param nodf_summary Output of `index_at_ages_summary` when index = "NODF".
#' @param q_sampled Output of `index_at_ages_samples` when index = "Q".
#' @param q_summary Output of `index_at_ages_summary` when index = "Q".
#' @param col_sampled Color used to represent values from sampled networks.
#' @param col_summary Color used to represent values from summary networks.
#'
#' @return A plot of z-scores over time. Violins show the posterior distribution of z-scores of sampled networks, with dots and lines showing the mean values. Z-scores of summary networks and extant network are shown as dots and line. Use different colors to differentiate values from sampled and summary networks.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # read data that comes with the package
#' data_path <- system.file("extdata", package = "evolnets")
#' tree <- read_tree_from_revbayes(paste0(data_path,"/tree_pieridae.tre"))
#' host_tree <- ape::read.tree(paste0(data_path,"/host_tree_pieridae.phy"))
#' history <- read_history(paste0(data_path,"/history_thin_pieridae.txt"), burnin = 0)
#'
#' # calculate posterior probabilities at ages
#' ages <- c(60, 50, 40, 0)
#' at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#'
#' # summary networks
#' summary_networks <- get_summary_networks(at_ages, threshold = 0.5, weighted = TRUE)
#' Nz_sum <- index_at_ages_summary(summary_networks, index = "NODF", nnull = 10)
#'
#' # sampled networks
#' sampled_networks <- get_sampled_networks(at_ages)
#' Nz_sam <- index_at_ages_samples(sampled_networks, index = "NODF", nnull = 10)
#'
#' # plot
#' plot_index_at_ages(nodf_sampled = Nz_sam, nodf_summary = Nz_sum)
plot_index_at_ages <- function(nodf_sampled, q_sampled = NULL, nodf_summary = NULL, q_summary = NULL, col_sampled = "#3B9AB2", col_summary = "#E67D00"){

  # Input checking
  # arguments must be data frames with columns `age`, `z` and `p` (only for samples)

  # Nestedness
  if(is.null(nodf_summary)) {
    max_z_n <- max(nodf_sampled$z, na.rm = TRUE)
  } else{
    max_z_n <- max(c(nodf_sampled$z, nodf_summary$z), na.rm = TRUE)
  }

  ppN <- nodf_sampled %>%
    dplyr::group_by(.data$age) %>%
    dplyr::filter(.data$p <= 0.05) %>%
    dplyr::summarise(pp = round(n()/max(.data$sample), digits = 2)) %>%
    dplyr::mutate(y = floor(max_z_n) + 3) # just for placement in the plot

  plotN <- ggplot2::ggplot(nodf_sampled) +
    ggplot2::geom_violin(ggplot2::aes(.data$age, .data$z, group = .data$age, fill = "Sampled"), col = NA, alpha = 0.5) +
    ggplot2::stat_summary(ggplot2::aes(.data$age, .data$z, col = "Sampled"), fun = "mean", geom = "line") +
    ggplot2::stat_summary(ggplot2::aes(.data$age, .data$z, col = "Sampled"), fun = "mean", geom = "point") +
    ggplot2::geom_text(ggplot2::aes(.data$age, .data$y, label = .data$pp), data = ppN) +
    ggplot2::scale_color_manual(name = "Network type",
                     breaks = c("Sampled", "Summary"),
                     values = c("Sampled" = col_sampled, "Summary" = col_summary)) +
    ggplot2::scale_fill_manual(breaks = c("Sampled"),
                     values = c("Sampled" = col_sampled)) +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma", fill = "") +
    ggplot2::theme_bw()

  if(!is.null(nodf_summary)) {
    plotN <- plotN +
      ggplot2::geom_point(ggplot2::aes(.data$age,.data$z, col = "Summary"), data = nodf_summary) +
      ggplot2::geom_line(ggplot2::aes(.data$age,.data$z, col = "Summary"), data = nodf_summary)

  }

  # Modularity
  if(!is.null(q_sampled)) {

    if(is.null(q_summary)) {
      max_z_q <- max(q_sampled, na.rm = TRUE)
    } else{
      max_z_q <- max(c(q_summary$z, q_sampled$z), na.rm = TRUE)
    }

    ppQ <- q_sampled %>%
      dplyr::group_by(.data$age) %>%
      dplyr::filter(.data$p <= 0.05) %>%
      dplyr::summarise(pp = round(n()/max(.data$sample), digits = 2))  %>%
      dplyr::mutate(y = floor(max_z_q) + 3)

    plotQ <- ggplot2::ggplot(q_sampled) +
      ggplot2::geom_violin(ggplot2::aes(.data$age, .data$z, group = .data$age, fill = "Sampled"), col = NA, alpha = 0.5) +
      ggplot2::stat_summary(ggplot2::aes(.data$age, .data$z, col = "Sampled"), fun = "mean", geom = "line") +
      ggplot2::stat_summary(ggplot2::aes(.data$age, .data$z, col = "Sampled"), fun = "mean", geom = "point") +
      ggplot2::geom_text(ggplot2::aes(.data$age, .data$y, label = .data$pp), data = ppQ) +
      ggplot2::scale_color_manual(name = "Network type",
                     breaks = c("Sampled", "Summary"),
                     values = c("Sampled" = col_sampled, "Summary" = col_summary)) +
      ggplot2::scale_fill_manual(breaks = c("Sampled"),
                     values = c("Sampled" = col_sampled)) +
      ggplot2::scale_x_reverse() +
      ggplot2::labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma", fill = "") +
      ggplot2::theme_bw()

    if(!is.null(q_summary)) {
      plotQ <- plotQ +
        ggplot2::geom_point(ggplot2::aes(.data$age,.data$z, col = "Summary"), data = q_summary) +
        ggplot2::geom_line(ggplot2::aes(.data$age,.data$z, col = "Summary"), data = q_summary)
    }

    return(patchwork::wrap_plots(plotN, plotQ, nrow = 2))
  } else {
    return(plotN)
  }
}


#' Network structure indices across ancestral summary networks and extant network
#'
#' Calculate z-scores for nestedness (NODF) or modularity (Q) for each ancestral network at time points in
#' the past and the extant network, based on null networks where all interactions have the same probability. By calculating
#' z-scores, we can compare ancestral networks at different ages.
#'
#' @param summary_networks List of ancestral networks at given ages and extant network. Usually from
#'   the output of `get_summary_networks`, see example.
#' @param index Index to be calculated for each ancestral network. "NODF" to calculate nestedness or
#'   "Q" to calculate modularity. For weighted networks, "weighted NODF" will be calculated.
#' @param ages Vector of ages (time points in the past) of ancestral networks. By default, uses all
#'   ages present in `summary_networks`.
#' @param nnull Number of null networks to generate to calculate the z-score. Default is 100.
#' @param use_future Parallelize with package `future`? Logical. Only applicable to modularity calculation.
#'
#' @return A data.frame of z-scores and p-values across networks.
#' @importFrom magrittr %>%
#' @importFrom methods is slot
#' @export
#'
#' @examples
#' # read data that comes with the package
#' data_path <- system.file("extdata", package = "evolnets")
#' tree <- read_tree_from_revbayes(paste0(data_path,"/tree_pieridae.tre"))
#' host_tree <- ape::read.tree(paste0(data_path,"/host_tree_pieridae.phy"))
#' history <- read_history(paste0(data_path,"/history_thin_pieridae.txt"))
#'
#' # get ancestral networks at ages in the past
#' ages <- c(60, 50, 40, 0)
#' at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#' summary_networks <- get_summary_networks(at_ages, threshold = 0.5, weighted = TRUE)
#'
#' # calculate nestedness of ancestral and extant networks
#' Nz <- index_at_ages_summary(summary_networks, index = "NODF")
#'
#' # calculate modularity of ancestral and extant networks with parallelization (slower)
#' # Qz <- index_at_ages_summary(summary_networks, index = "Q", use_future = TRUE)
index_at_ages_summary <- function(summary_networks, index, ages = NULL, nnull = 100, use_future= FALSE){

  # Input checking
  if (!is.list(summary_networks) || !all(vapply(summary_networks, inherits, TRUE, 'matrix'))) {
    stop('`summary_networks` should be a list of matrices, usually generated by `get_summary_networks`.')
  }
  index <- match.arg(index, c('NODF', 'Q'), several.ok = FALSE)
  if (!is.null(ages) & !is.numeric(ages)) stop('`ages` should be numeric.')
  if (!is.numeric(nnull)) stop('`nnull` should be numeric.')

  if (is.null(ages)) ages <- as.numeric(names(summary_networks))

  # Keep only specified ages
  summary_networks <- summary_networks[as.character(ages)]

  # Calculating indices
  if (index == "NODF") {

    NODF <- tibble::tibble()

    for(a in seq_along(ages)){
      network <- summary_networks[[a]]

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(network) < 2 || nrow(network) < 2 || is.vector(network)) {
        warning(paste0("Skipping network at age ",ages[a]," because it has less than 2 hosts or symbionts"))
      } else {
        NODF_age <- get_z_nodf(network, nnull) %>%
          dplyr::mutate(age = ages[a])
        NODF <- dplyr::bind_rows(NODF, NODF_age)
      }
    }
    ret <- data.frame(NODF)
  }

  if (index == "Q") {

    Q <- tibble::tibble()

    for(a in seq_along(ages)){
      network <- summary_networks[[a]]

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(network) < 2 || nrow(network) < 2 || is.vector(network)) {
        warning(paste0("Skipping network at age ",ages[a],"because it has less than 2 hosts or symbionts"))
      } else {
        Q_age <- get_z_q(network, nnull, use_future = use_future)  %>%
          dplyr::mutate(age = ages[a])
        Q <- dplyr::bind_rows(Q, Q_age)
      }
    }
    ret <- data.frame(Q)
  }

  return(ret)
}


# calculate z-score for nestedness of a given extant or ancestral network
get_z_nodf <- function(network, nnull = 100){

  values <- mapply(unique, network) %>%
    unique() %>%
    sort() %>%
    as.numeric()

  # if network is weighted or has both 1s and 2s, use weighted NODF
  if(identical(values, c(0,1))) {

    # calculate NODf for observed network
    Nobs <- bipartite::networklevel(network, index="NODF")

    # generate null networks
    null_model <- vegan::nullmodel(network, "r00")
    null_nets <- stats::simulate(null_model, nsim=nnull)

    # calculate NODF

    Nnull <- tibble::tibble()

    for(j in 1:nnull){
      Nrandom <- bipartite::networklevel(null_nets[,,j], index="NODF")
      Nnull <- dplyr::bind_rows(Nnull, tibble(null_idx = j, NODF = Nrandom))
    }

  } else{
    # calculate NODF for observed network
    Nobs <- bipartite::networklevel(network, index="weighted NODF")

    # generate null networks
    count <- round(network*100)            # transform the probabilities into counts
    # to use the null model 'r00_both'
    null_model <- vegan::nullmodel(count, "r00_both")
    null_nets <- stats::simulate(null_model, nsim=nnull)

    # calculate NODF for null networks
    Nnull <- tibble::tibble()

    for(j in 1:nnull){
      Nrandom <- bipartite::networklevel(null_nets[,,j], index="weighted NODF")
      Nnull <- dplyr::bind_rows(Nnull, tibble(null_idx = j, NODF = Nrandom))
    }

  }

  Nz <- Nnull %>%
    summarize(mean = mean(.data$NODF),
              sd = stats::sd(.data$NODF)) %>%
    mutate(z = (Nobs - .data$mean)/.data$sd,
           N_obs = Nobs) %>%
    dplyr::select(.data$N_obs, .data$mean, .data$sd, .data$z)

  return(data.frame(Nz))

}


# calculate z-score for modularity of a given extant or ancestral network
get_z_q <- function(network, nnull = 100, use_future=FALSE){

  values <- mapply(unique, network) %>%
    unique() %>%
    sort() %>%
    as.numeric()

  # generate null networks
  # if network is weighted or has both 1s and 2s, use r00_both
  if(identical(values, c(0,1))) {
    null_model <- vegan::nullmodel(network, "r00")
    null_nets <- stats::simulate(null_model, nsim=nnull)
  } else{
    count <- round(network*100)            # transform the probabilities into counts
                                           # to use the null model 'r00_both'
    null_model <- vegan::nullmodel(count, "r00_both")
    null_nets <- stats::simulate(null_model, nsim=nnull)
  }

  # calculate Q
  Qnull <- tibble::tibble()

  # use future
  if (use_future) {
    Qrandom <- do.call(c, future.apply::future_lapply(1:nnull, function(x) { mycomputeModules(null_nets[,,x])@likelihood }, future.seed=T))
    Qnull <- tibble::tibble(null_idx = 1:nnull, Q = Qrandom)

  }
  # use base/serial
  else {
    for(j in 1:nnull){
      Qrandom <- mycomputeModules(null_nets[,,j])@likelihood
      Qnull <- dplyr::bind_rows(Qnull, tibble(null_idx = j, Q = Qrandom))
    }
  }

  # calculate Q for observed network
  Qobs <- mycomputeModules(network)@likelihood

  Qz <- Qnull %>%
    summarize(mean = mean(.data$Q),
              sd = stats::sd(.data$Q)) %>%
    mutate(z = (Qobs - .data$mean)/.data$sd,
           Q_obs = Qobs) %>%
    dplyr::select(.data$Q_obs, .data$mean, .data$sd, .data$z)

  return(Qz)

}


#' Posterior distribution of network structure indices across ancestral networks
#'
#' Calculate z-scores for nestedness (NODF) or modularity (Q) for each MCMC sample at time points in
#' the past based on null networks where all interactions have the same probability. By calculating
#' z-scores, we can compare ancestral networks at different ages.
#'
#' @param sampled_networks List of ancestral networks sampled across MCMC at given ages. Usually from
#'   the output of `get_sampled_networks`, see example.
#' @param index Index to be calculated for each ancestral network. "NODF" to calculate nestedness or
#'   "Q" to calculate modularity.
#' @param ages Vector of ages (time points in the past) of ancestral networks. By default, uses all
#'   ages present in `sampled_networks`.
#' @param nnull Number of null networks to generate to calculate the z-score. Default is 100.
#' @param use_future Parallelize with package `future`? Logical. Only applicable to modularity calculation.
#'
#' @return A data.frame of z-scores and p-values across samples and ages.
#' @importFrom magrittr %>%
#' @importFrom methods is slot
#' @export
#'
#' @examples
#' # read data that comes with the package
#' data_path <- system.file("extdata", package = "evolnets")
#' tree <- read_tree_from_revbayes(paste0(data_path,"/tree_pieridae.tre"))
#' host_tree <- ape::read.tree(paste0(data_path,"/host_tree_pieridae.phy"))
#' history <- read_history(paste0(data_path,"/history_thin_pieridae.txt"), burnin = 0)
#'
#' # get sampled networks at ages in the past
#' ages <- c(60,50,40,0)
#' samples_at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#' sampled_networks <- get_sampled_networks(samples_at_ages)
#'
#' # calculate posterior distribution of nestedness
#' Nz <- index_at_ages_samples(sampled_networks, index = "NODF")
#'
#' #  calculate posterior distribution of modularity with parallelization (slower)
#' # Qz <- index_at_ages_samples(sampled_networks, index = "Q", use_future = TRUE)
index_at_ages_samples <- function(sampled_networks, index, ages = NULL, nnull = 100, use_future = FALSE){

  # Input checking
  if (!is.list(sampled_networks)) {
    stop('`sampled_networks` should be a list, usually generated by `get_sampled_networks`.')
  }
  if (!(all(vapply(sampled_networks, is.array, TRUE)))) {
    stop('All list entries of `sampled_networks` should be arrays.')
  }
  index <- match.arg(index, c('NODF', 'Q'), several.ok = FALSE)
  if (!is.null(ages) & !is.numeric(ages)) stop('`ages` should be numeric.')
  if (!is.numeric(nnull)) stop('`nnull` should be numeric.')

  if (is.null(ages)) ages <- as.numeric(names(sampled_networks))

  # Keep only specified ages
  sampled_networks <- sampled_networks[as.character(ages)]

  # Calculating indices
  if (index == "NODF") {
    # find unique values, to check for three-state model.
    unique_vals <- Reduce(union, lapply(sampled_networks, function(x) unique(c(x))))
    three_state <- identical(sort(as.numeric(unique_vals)), 0:2)

    NODF_null <- NODF_samples_null(sampled_networks, ages, nnull, weighted = three_state)
    NODF_samples <- NODF_samples_at_ages(sampled_networks, ages, weighted = three_state)

    NODF_pvals <- NODF_null %>%
      dplyr::filter(!is.na(.data$NODFnull)) %>%
      dplyr::left_join(NODF_samples, by = c("age", "sample")) %>%
      dplyr::group_by(.data$age, .data$sample) %>%
      dplyr::summarise(p = sum(.data$NODFnull >= .data$obs_NODF) / nnull, .groups = 'drop')

    NODF_zsamples <- NODF_null %>%
      dplyr::group_by(.data$age, .data$sample) %>%
      dplyr::summarize(
        mean_NODF = mean(.data$NODFnull),
        sd_NODF = stats::sd(.data$NODFnull),
        .groups = 'drop'
      ) %>%
      dplyr::left_join(NODF_samples, by = c("age", "sample")) %>%
      dplyr::mutate(z = (.data$obs_NODF - .data$mean_NODF) / .data$sd_NODF) %>%
      dplyr::left_join(NODF_pvals, by = c("age", "sample"))

    ret <- as.data.frame(NODF_zsamples)
  }

  if (index == "Q") {

    Q_null <- Q_samples_null(sampled_networks, ages, nnull, use_future = use_future)
    Q_samples <- Q_samples_at_ages(sampled_networks, ages)

    Q_pvals <- Q_null %>%
      dplyr::filter(!is.na(.data$Qnull)) %>%
      dplyr::left_join(Q_samples, by = c("age", "sample")) %>%
      dplyr::group_by(.data$age, .data$sample) %>%
      dplyr::summarise(p = sum(.data$Qnull >= .data$obs_Q) / nnull, .groups = 'drop')

    Q_zsamples <- Q_null %>%
      dplyr::group_by(.data$age, .data$sample) %>%
      dplyr::summarize(
        mean_Q = mean(.data$Qnull),
        sd_Q = stats::sd(.data$Qnull),
        .groups = 'drop'
      ) %>%
      dplyr::left_join(Q_samples, by = c("age", "sample")) %>%
      dplyr::mutate(z = (.data$obs_Q - .data$mean_Q) / .data$sd_Q) %>%
      dplyr::left_join(Q_pvals, by = c("age", "sample"))

    ret <- as.data.frame(Q_zsamples)
  }

  ret <- ret[!is.na(ret[,1]),]
  return(ret)
}


# Simulate null networks and calculate NODF
NODF_samples_null <- function(sampled_networks, ages, nnull, weighted = FALSE){
  index <- ifelse(weighted, 'weighted NODF', 'NODF')
  null_type <- ifelse(weighted, 'r00_both', 'r00')

  nsamp <- dim(sampled_networks[[1]])[1]

  NODF_null <- data.frame(
    age = rep.int(NA, nnull * nsamp * length(sampled_networks)), sample = NA, sim = NA, NODFnull = NA
  )

  for (a in seq_along(sampled_networks)) {
    for (i in seq_len(dim(sampled_networks[[a]])[1])) {
      net <- sampled_networks[[a]][i, , ]
      net <- net[rowSums(net) != 0, ]
      net <- net[, colSums(net) != 0]

      nullm <- vegan::nullmodel(net, null_type)
      sim <- stats::simulate(nullm, nsim = nnull)

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(net) < 2 || nrow(net) < 2 || is.vector(net)) {
        warning(paste0("Skipping network at age ",ages[a]," because it has less than 2 hosts or symbionts"))
      } else {
        for (j in seq_len(nnull)) {
        Nrandom <- bipartite::networklevel(sim[, , j], index = index)
        NODF_null[(a - 1) * nsamp * nnull + (i - 1) * nnull + j, ] <- c(ages[a], i, j, Nrandom)
        }
      }
    }
  }

  NODF_null

}


# Get NODF for each sampled network
NODF_samples_at_ages <- function(sampled_networks, ages, weighted = FALSE){
  index <- ifelse(weighted, 'weighted NODF', 'NODF')

  nsamp <- dim(sampled_networks[[1]])[1]
  NODF_samples <- data.frame(
    age = rep.int(NA, nsamp * length(sampled_networks)), sample = NA, obs_NODF = NA
  )

  for (a in seq_along(sampled_networks)) {
    for (i in seq_len(dim(sampled_networks[[a]])[1])) {
      net <- sampled_networks[[a]][i, , ]
      net <- net[rowSums(net) != 0, ]
      net <- net[, colSums(net) != 0]

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(net) < 2 || nrow(net) < 2 || is.vector(net)) {
        warning(paste0("Skipping network at age ",ages[a]," because it has less than 2 hosts or symbionts"))
      } else {
        nodf <- bipartite::networklevel(net, index = index)
        NODF_samples[(a - 1) * nsamp + i, ] <- c(ages[a], i, nodf)
      }
    }
  }

  NODF_samples

}

# Simulate null networks and calculate Q
Q_samples_null <- function(sampled_networks, ages, nnull, use_future = FALSE){

  nsamp <- dim(sampled_networks[[1]])[1]
  Q_null <- data.frame(
    age = rep.int(NA, nnull * nsamp * length(sampled_networks)), sample = NA, sim = NA, Qnull = NA
  )

  for (a in seq_along(sampled_networks)) {
    for (i in seq_len(dim(sampled_networks[[a]])[1])) {
      net <- sampled_networks[[a]][i, , ]
      net <- net[rowSums(net) != 0, ]
      net <- net[, colSums(net) != 0]

      nullm <- vegan::nullmodel(net, "r00")
      sim <- stats::simulate(nullm, nsim = nnull)

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(net) < 2 || nrow(net) < 2 || is.vector(net)) {
        warning(paste0("Skipping network at age ",ages[a]," because it has less than 2 hosts or symbionts"))
      } else {
        if (use_future) {
          Qrandom <- do.call(c, future.apply::future_lapply(1:nnull, function(x) { mycomputeModules(sim[,,x])@likelihood }, future.seed=T))
          for (j in 1:nnull) {
            Q_null[(a - 1) * nsamp * nnull + (i - 1) * nnull + j, ] <- c(ages[a], i, j, Qrandom[j])
          }
        } else {
          for (j in 1:nnull) {
            mod <- mycomputeModules(sim[, , j])
            Qrandom <- mod@likelihood
            Q_null[(a - 1) * nsamp * nnull + (i - 1) * nnull + j, ] <- c(ages[a], i, j, Qrandom)
          }
        }
      }
    }
  }

  Q_null

}


# Get Q for each sampled network
Q_samples_at_ages <- function(sampled_networks, ages){

  nsamp <- dim(sampled_networks[[1]])[1]
  Q_samples <- data.frame(
    age = rep.int(NA, nsamp * length(sampled_networks)), sample = NA, obs_Q = NA
  )

  for (a in seq_along(sampled_networks)) {
    for (i in seq_len(dim(sampled_networks[[a]])[1])) {
      net <- sampled_networks[[a]][i,,]
      net <- net[rowSums(net) != 0, ]
      net <- net[, colSums(net) != 0]

      # Skip networks with less than 2 hosts or 2 symbionts
      if(ncol(net) < 2 || nrow(net) < 2 || is.vector(net)) {
        warning(paste0("Skipping network at age ",ages[a]," because it has less than 2 hosts or symbionts"))
      } else {
        mod <- mycomputeModules(net)
        Q <- mod@likelihood
        Q_samples[(a - 1) * nsamp + i, ] <- c(ages[a], i, Q)
      }
    }
  }

  Q_samples

}
