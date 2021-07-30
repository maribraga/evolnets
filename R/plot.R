#' Plot a network with modules modules as an adjacency matrix
#'
#' @param net An adjacency matrix for a bipartite network. Parasites should be the rows, hosts
#'   should be columns. If all values are 0 or 1 an unweighted network is represented, otherwise
#'   a weighted network is assumed.
#' @param modules A `moduleWeb` object defining the models in the network. If left `NULL` the
#'   modules are automatically calculated.
#' @param parasite_order A character vector giving the order the parasite should be listed in.
#'   Should contain each parasite only once, and include the row names of `net`.
#' @param host_order A character vector giving the order the hosts should be listed in. Should
#'   contain each host only once, and include the column names of `net`.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # The slow portion of this function is the calculation of the modules.
#'   plot_module_matrix(extant_net)
#'
#'   # Change our network to a weighted one:
#'   extant_net_weighted <- extant_net
#'   extant_net_weighted[extant_net == 1] <- runif(sum(extant_net))
#'   plot_module_matrix(extant_net_weighted)
#' }
plot_module_matrix <- function(net, modules = NULL, parasite_order = NULL, host_order = NULL) {
  # Check inputs.
  if (!is.matrix(net)) stop('`net` should be a matrix.')
  if (!is.null(modules) && !inherits(modules, 'moduleWeb')) {
    stop('`modules` should be of class `moduleWeb`.')
  }
  if (!is.null(parasite_order) && (
    length(setdiff(rownames(net), parasite_order)) != 0 || !is.vector(parasite_order) ||
    !is.character(parasite_order) || any(duplicated(parasite_order))
  )) {
    stop('`parasite_order` should be a character vector without duplicates.')
  }
  if (!is.null(host_order) && (
    length(setdiff(colnames(net), host_order)) != 0 || !is.vector(host_order) ||
    !is.character(host_order) || any(duplicated(host_order))
  )) {
    stop('`host_order` should be a character vector without duplicates.')
  }

  # If no modules are given, calculate them
  if (is.null(modules)) {
    modules <- mycomputeModules(net)
  }
  # Take the modules, and put the host modules in a data.frame
  mod_list <- listModuleInformation(modules)[[2]]
  host_mods <- lapply(mod_list, function(x) data.frame(host = x[[2]]))
  host_mods <- dplyr::bind_rows(host_mods, .id = 'host_module')
  para_mods <- lapply(mod_list, function(x) data.frame(parasite = x[[1]]))
  para_mods <- dplyr::bind_rows(para_mods, .id = 'parasite_module')

  # Take the extant network (as an adjacency matrix), and create a plottable data.frame
  net_df <- as.data.frame(net)
  net_df$parasite <- rownames(net_df)
  net_df <- tidyr::pivot_longer(net_df, -.data$parasite, names_to = 'host', values_to = 'weight')
  net_df <- net_df[net_df$weight != 0, ]

  # Join the extant network with the module info
  module_mat <- dplyr::left_join(net_df, host_mods, by = 'host')
  module_mat <- dplyr::left_join(module_mat, para_mods, by = 'parasite')
  module_mat$module <- ifelse(
    module_mat$host_module == module_mat$parasite_module,
    module_mat$host_module,
    NA
  )
  # Then join with the order of tips from the trees to ensure correct alignment
  if (!is.null(parasite_order)) {
    module_mat$parasite <- factor(module_mat$parasite, levels = parasite_order)
  }
  if (!is.null(host_order)) {
    module_mat$host <- factor(module_mat$host, levels = host_order)
  }
  if (all(module_mat$weight %in% c(0L, 1L))) {
    p <- ggplot2::ggplot(module_mat, ggplot2::aes_(~host, ~parasite, fill = ~module))
  } else {
    p <- ggplot2::ggplot(module_mat, ggplot2::aes_(~host, ~parasite, fill = ~module, alpha = ~weight))
  }

  p +
    ggplot2::geom_hline(yintercept = 0.5 + 0:nrow(module_mat), color = 'grey80', size = 0.1) +
    ggplot2::geom_vline(xintercept = 0.5 + 0:nrow(module_mat), color = 'grey80', size = 0.1) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0.5)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0.5)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5),
      axis.ticks = ggplot2::element_blank(),
      legend.position = 'top'
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = 'Module'
    )
}

#' Plot ancestral states on the phylogeny.
#'
#' @param tree The phylogeny, a `phylo` object.
#' @param samples_at_nodes A list of length 2, output from `posterior_at_nodes()`.
#' @param modules A `moduleWeb` object defining the models in the network.
#' @param threshold The posterior probability above which the ancestral states should be shown.
#'   Defaults to 90% (`0.9`). Numeric vector of length 1.
#' @param point_size How large the ancestral state points should be, default at 3. Play with this
#'   and `dodge_width` to get a pleasing result. Numeric vector of length 1.
#' @param dodge_width How far the points should be pushed apart, when there is multiple states at
#'   a single node, default at 0.025. Play with this and `point_size` to get a pleasing result.
#'   Numeric vector of length 1.
#' @param legend Whether to display a legend for the colors. Logical vector of length 1.
#' @param colors Override the default colors. Should be a character vector with as many color values
#'  as there are modules.
#'
#' The ancestral states are automatically colored by module. To change what colors are used, you
#' can add color scales to the resulting `ggplot`, e.g. `ggplot2::scale_color_manual()`.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   san <- posterior_at_nodes(history, 66 + 1:65, host_tree)
#'   mods <- mycomputeModules(extant_net)
#'   plot_ancestral_states(tree, san, mods)
#'   # Manual colors
#'   plot_ancestral_states(tree, san, mods, colors = rainbow(20))
#' }
plot_ancestral_states <- function(
  tree, samples_at_nodes, modules,
  threshold = 0.9, point_size = 3, dodge_width = 0.025, legend = TRUE, colors = NULL
) {
  if (!requireNamespace('ggtree')) {
    stop('Please install the ggtree package to use this function. Use `BiocManager::install("ggtree")`')
  }
  # Check inputs
  if (!inherits(tree, 'phylo')) stop('`tree` should be of class `phylo`.')
  if (!is.list(samples_at_nodes)) stop('`samples_at_nodes` should be a list.')
  if (!all(rownames(samples_at_nodes[[2]]) %in% tree$node.label)) {
    stop('All nodes in `samples_at_nodes` should occur in the tree as internal nodes, check `tree$node.label`')
  }

  # Take the modules, and put the host modules in a data.frame
  mod_list <- listModuleInformation(modules)[[2]]
  host_mods <- lapply(mod_list, function(x) data.frame(host = x[[2]]))
  host_mods <- dplyr::bind_rows(host_mods, .id = 'module')

  # Get the ancestral states, reformat to plot on the parasite tree
  node_df <- as.data.frame(samples_at_nodes[[2]])
  node_df$node <- rownames(node_df)
  node_df <- tidyr::pivot_longer(node_df, -.data$node, names_to = 'host')
  node_df <- node_df[node_df$value >= threshold, ]
  node_df$value <- NULL
  node_df <- dplyr::left_join(node_df, host_mods, by = "host")

  # Match the ancestral states with the right plotting coordinates
  # Little helper function to manually dodge the points.
  offset_x <- function(x, width) {
    stopifnot(all(x == x[1]))
    hw <- floor(length(x) / 2)
    if (length(x) %% 2 == 0) {
      out <- x + width *
        c(rev(seq.int(-0.5, by = -1, length.out = hw)), seq.int(0.5, by = 1, length.out = hw))
    } else {
      out <- x + width *
        c(rev(seq.int(-1, by = -1, length.out = hw)), 0, seq.int(1, by = 1, length.out = hw))
    }
    return(out)
  }
  # Do we need a legend?
  if (legend) guide <- 'legend' else guide <- 'none'
  # Set up our color scale
  mods <- seq_along(mod_list)
  if (is.null(colors)) {
    color_scale <- ggplot2::scale_color_discrete(
      limits = factor(mods, levels = mods), guide = guide
    )
  } else {
    color_scale <- ggplot2::scale_color_manual(
      values = colors, limits = factor(mods, levels = mods), guide = guide
    )
  }

  # Make the parasite tree
  p <- ggplot2::ggplot(tree) +
    ggtree::geom_tree() +
    ggplot2::scale_x_continuous(name = NULL, labels = abs, expand = ggplot2::expansion(c(0.05, 0))) +
    ggplot2::scale_y_continuous(expand = c(0, 0.5)) +
    color_scale +
    ggtree::theme_tree2()
  # Flip the time axis the right way around
  p <- ggtree::revts(p)

  # Extract the node coordinates, so we can easily plot our own node information
  coords <- p$data[, c('x', 'y', 'label', 'isTip')]
  coords <- coords[order(coords$y), ]

  x_range <- diff(range(coords$x))
  node_df2 <- dplyr::left_join(node_df, coords, by = c('node' = 'label')) %>%
    dplyr::arrange(.data$node, .data$module) %>%
    dplyr::group_by(.data$node) %>%
    dplyr::mutate(
      x = offset_x(.data$x, width = x_range * dodge_width)
    )
  p <- p + ggplot2::geom_point(
      ggplot2::aes_(~x, ~y, color = ~module),
      node_df2, shape = 15, size = point_size
    )

  return(p)
}

#' Plot a network with modules modules as an adjacency matrix, with aligned phylogenies.
#'
#' @param net An adjacency matrix for a bipartite network. This should be the extant network.
#'   Parasites should be the rows, hosts should be columns. If all values are 0 or 1 an unweighted
#'   network is represented, otherwise a weighted network is assumed.
#' @param tree The phylogeny belonging the parasites or herbivores. Object of class `phylo`.
#' @param host_tree The phylogeny belonging to the hosts. Object of class `phylo`.
#'
#' See the examples on how to change the color scale.
#'
#' @return An assembly of plots, of class `patchwork`.
#' @export
#'
#' @inheritParams plot_ancestral_states
#' @inheritParams plot_module_matrix
#'
#' @examples
#' \dontrun{
#'   san <- posterior_at_nodes(history, 66 + 1:65, host_tree)
#'   plot_module_matrix2(extant_net, san, tree, host_tree)
#'   # manual_colors
#'   plot_module_matrix2(extant_net, san, tree, host_tree, colors = rainbow(20))
#' }
plot_module_matrix2 <- function(
  net, samples_at_nodes, tree, host_tree,
  modules = NULL, threshold = 0.9, point_size = 3, dodge_width = 0.025, colors = NULL
) {
  # Check inputs
  if (!is.matrix(net)) stop('`net` should be a matrix.')
  if (!inherits(host_tree, 'phylo')) stop('`host_tree` should be of class `phylo`.')
  # Other inputs will be checked by downstream functions.

  # Get module information for the net, if not provided
  if (is.null(modules)) {
    modules <- mycomputeModules(net)
  }

  # Make the parasite tree plot
  parasite_plot <- plot_ancestral_states(
    tree, samples_at_nodes, modules,
    threshold = threshold, point_size = point_size, dodge_width = dodge_width, legend = FALSE,
    colors = colors
  )
  parasite_coords <- parasite_plot$data[parasite_plot$data$isTip, c('x', 'y', 'label', 'isTip')]
  parasite_coords <- parasite_coords[order(parasite_coords$y), ]

  # Make the host tree plot
  mods <- seq_along(listModuleInformation(modules)[[2]])
  host_plot <- ggplot2::ggplot(host_tree) +
    ggtree::geom_tree() +
    ggplot2::coord_flip() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0.5)) +
    ggtree::theme_tree()
  host_coords <- host_plot$data[host_plot$data$isTip, c('y', 'label')]
  host_coords <- host_coords[order(host_coords$y), ]

  # Make the matrix
  module_plot <- plot_module_matrix(net, modules, parasite_coords$label, host_coords$label)
  if (is.null(colors)) {
    module_plot <- module_plot + ggplot2::scale_fill_discrete(limits = factor(mods, levels = mods))
  } else {
    module_plot <- module_plot +
      ggplot2::scale_fill_manual(values = colors, limits = factor(mods, levels = mods))
  }

  patchwork::wrap_plots(
    parasite_plot + ggplot2::labs(tag = 'A'),
    module_plot + ggplot2::labs(tag = 'B'),
    patchwork::guide_area(),
    host_plot + ggplot2::labs(tag = 'C'),
    ncol = 2, widths = c(1, 1), heights = c(3, 1), guides = 'collect'
  )
}
