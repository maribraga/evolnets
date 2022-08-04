#' Plot a network with modules as an adjacency matrix
#'
#' @param net An adjacency matrix for a bipartite network. Parasites should be
#'   the rows, hosts should be columns. If all values are 0 or 1 an unweighted
#'   network is represented, otherwise a weighted network is assumed.
#' @param modules A `moduleWeb` or a `data.frame` object defining the models in
#'   the network. If left `NULL` the modules are automatically calculated. If a
#'   `data.frame` is passed, it must contain three columns: $name with taxon
#'   names, $module with the module the taxon was assigned to, and $type which
#'   defines if the taxon is a "host" or a "symbiont".
#' @param parasite_order A character vector giving the order the parasite should
#'   be listed in. Should contain each parasite only once, and include the row
#'   names of `net`.
#' @param host_order A character vector giving the order the hosts should be
#'   listed in. Should contain each host only once, and include the column names
#'   of `net`.
#' @param module_order A character vector giving the order that modules should
#'   be plotted. Should contain each module only once.
#' @param state_alpha A numeric vector of length 2. Gives the alpha
#'   (transparency) values for the interaction type in the three-state model
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # The slow portion of this function is the calculation of the modules.
#'   plot_extant_matrix(extant_net)
#'
#'   # Change our network to a weighted one:
#'   extant_net_weighted <- extant_net
#'   extant_net_weighted[extant_net == 1] <- runif(sum(extant_net))
#'   plot_extant_matrix(extant_net_weighted)
#' }
plot_extant_matrix <- function(
    net, modules = NULL, module_order = NULL, parasite_order = NULL,
    host_order = NULL, state_alpha = c(0.5, 1)
) {
  # Check inputs.
  if (!is.matrix(net)) stop('`net` should be a matrix.')
  if (!is.null(modules) && (
    !inherits(modules, 'moduleWeb') && !inherits(modules, 'data.frame')
  )) {
    stop('`modules` should be of class `moduleWeb` or `data.frame`.')
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
  if (!is.numeric(state_alpha) || length(state_alpha) != 2) {
    stop('`state_alpha` should be a numeric vector of length 2.')
  }

  # If no modules are given, calculate them
  if (is.null(modules)) {
    modules <- mycomputeModules(net)
  }

  if (inherits(modules, 'moduleWeb')) {
    # Take the modules, and put the host and symbiont modules in data.frames
    mod_list <- bipartite::listModuleInformation(modules)[[2]]
    host_mods <- lapply(mod_list, function(x) data.frame(host = x[[2]]))
    host_mods <- dplyr::bind_rows(host_mods, .id = 'host_module')
    host_mods$host_module <- as.numeric(host_mods$host_module)
    para_mods <- lapply(mod_list, function(x) data.frame(parasite = x[[1]]))
    para_mods <- dplyr::bind_rows(para_mods, .id = 'parasite_module')
    para_mods$parasite_module <- as.numeric(para_mods$parasite_module)
  } else {
    if (inherits(modules, 'data.frame') && !all(
      c('name', 'module', 'type') %in% names(modules)
    )) {
      stop('`modules` should be a `moduleWeb` or a data.frame with "name", "module" and "type" columns.')
    }
    host_mods <- modules %>%
      dplyr::filter(.data$type == "host") %>%
      dplyr::select(.data$name, .data$module) %>%
      dplyr::rename(host = .data$name, host_module = .data$module)
    para_mods <- modules %>%
      dplyr::filter(.data$type == "symbiont") %>%
      dplyr::select(.data$name, .data$module) %>%
      dplyr::rename(parasite = .data$name, parasite_module = .data$module)
  }

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
  # Set order of modules
  if (is.null(module_order)) {
    module_order <- sort(unique(module_mat$parasite_module))
  }
  module_mat$module <- factor(module_mat$module, levels = module_order)

  # Then join with the order of tips from the trees to ensure correct alignment
  if (is.null(parasite_order)) {
    parasite_order <- module_mat %>%
      dplyr::group_by(.data$parasite, .data$parasite_module) %>%
      dplyr::summarise(degree = dplyr::n()) %>%
      dplyr::mutate(parasite_module = factor(.data$parasite_module, levels = module_order)) %>%
      dplyr::arrange(.data$parasite_module,.data$degree) %>%
      dplyr::pull(.data$parasite)
  }
  if (is.null(host_order)) {
    host_order <- module_mat %>%
      dplyr::group_by(.data$host, .data$host_module) %>%
      dplyr::summarise(degree = dplyr::n()) %>%
      dplyr::mutate(host_module = factor(.data$host_module, levels = module_order)) %>%
      dplyr::arrange(.data$host_module, dplyr::desc(.data$degree)) %>%
      dplyr::pull(.data$host)
  }
  module_mat$parasite <- factor(module_mat$parasite, levels = parasite_order)
  module_mat$host <- factor(module_mat$host, levels = host_order)

  net_vals <- sort(unique(as.numeric(net)))

  if (identical(net_vals, c(0, 1, 2))) {
    p <- ggplot2::ggplot(module_mat, ggplot2::aes_(~host, ~parasite, fill = ~module, alpha = ~factor(weight))) +
      ggplot2::scale_alpha_ordinal(
        limits = factor(1:2), range = state_alpha, name = 'Interaction type',
        labels = c('1' = 'Potential', '2' = 'Actual'),
        guide = ggplot2::guide_legend(ncol = 1)
      )
  } else if (identical(net_vals, c(0, 1))) {
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
#' @param at_nodes A list of length 2, output from `posterior_at_nodes()`.
#' @param modules A `moduleWeb` or a `data.frame` object defining the modules in the network.
#' If a `data.frame` is passed, it must contain three columns:
#'   $name with taxon names,
#'   $module with the module the taxon was assigned to, and
#'   $type which defines if the taxon is a "host" or a "symbiont".
#' @param module_order A character vector giving the order that modules should be plotted. Should contain
#'   each module only once.
#' @param type One of `'states'` or `'repertoires'`. If `'states'`, will plot the presence of a
#'   state when its posterior probability is higher than `threshold`. If `'repertoires'`, will plot
#'   the same but for the given `repertoire`.
#' @param state Which state? Default is 2. For analyses using the 3-state model, choose `1` or `2`
#'    (where 1 is a potential host and 2 an actual host). Only used if `type` is `'states'`.
#' @param repertoire Either the `'realized'` repertoire which is defined as state 2, or the
##'  `'fundamental'` repertoire (default) which is defined as having any state (usually 1 or 2), and
##'  in the 3-state model includes both actual and potential hosts.
#' @param layout One of `'rectangular'`, `'slanted'`, `'fan'`, `'circular'`, `'radial'`,
#'   `'equal_angle'`, `'daylight'` or `'ape'`.
#' @param threshold The posterior probability above which the ancestral states should be shown.
#'   Defaults to 90% (`0.9`). Numeric vector of length 1.
#' @param point_size How large the ancestral state points should be, default at 3. Play with this
#'   and `dodge_width` to get a pleasing result. Numeric vector of length 1.
#' @param point_shape What point shape should be used for the ancestral states? When left `NULL`,
#'   a reasonable default will be chosen. Otherwise, a numeric vector of length 1.
#' @param dodge_width How far the points should be pushed apart, when there is multiple states at
#'   a single node, default at 0.025. Play with this and `point_size` to get a pleasing result.
#'   Numeric vector of length 1.
#' @param legend Whether to display a legend for the colors. Logical vector of length 1.
#' @param colors Override the default colors. Should be a character vector with as many color values
#'  as there are modules.
#' @param state_alpha A numeric vector of length 2. Gives the alpha
#'   (transparency) values for the interaction type in the three-state model
#' @param ladderize Logical. Whether to ladderize the tree. Default to FALSE.
#'
#' The ancestral states are automatically colored by module. To change what colors are used, you
#' can add color scales to the resulting `ggplot`, e.g. `ggplot2::scale_color_manual()`.
#'
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   at_nodes <- posterior_at_nodes(history, tree, host_tree, 66 + 1:65)
#'   mods <- mycomputeModules(extant_net)
#'   plot_ancestral_states(tree, at_nodes, mods)
#'   # Manual colors
#'   plot_ancestral_states(tree, at_nodes, mods, colors = rainbow(20))
#' }
plot_ancestral_states <- function(
    tree, at_nodes, modules, module_order = NULL, type = "states", state = 2,
    repertoire = 'fundamental', layout = "rectangular", threshold = 0.9, point_size = 3,
    point_shape = NULL, dodge_width = 0.025, legend = TRUE, colors = NULL,
    state_alpha = c(0.5, 1), ladderize = FALSE
) {
  if (!requireNamespace('ggtree')) {
    stop('Please install the ggtree package to use this function. Use `BiocManager::install("ggtree")`')
  }
  # Check inputs
  if (!inherits(tree, 'phylo')) stop('`tree` should be of class `phylo`.')
  if (!is.list(at_nodes)) stop('`at_nodes` should be a list.')
  if (!all(colnames(at_nodes[['samples']]) %in% tree$node.label)) {
    stop('All nodes in `at_nodes` should occur in the tree as internal nodes, check `tree$node.label`')
  }
  if (!(type %in% c('states', 'repertoires') & length(type) == 1)) {
    stop("`type` should be either 'states' or 'repertoires'.")
  }
  if (!all(as.character(state) %in% dimnames(at_nodes[['post_states']])[[3]])) {
    stop("`state` should be a vector and have valid states occuring in `at_nodes`")
  }
  if (!(repertoire %in% c('fundamental', 'realized') & length(repertoire) == 1)) {
    stop("`repertoire` should be either 'fundamental' or 'realized'.")
  }

  if (inherits(modules, 'moduleWeb')) {
    # Take the modules, and put the host modules in a data.frame
    mod_list <- bipartite::listModuleInformation(modules)[[2]]
    host_mods <- lapply(mod_list, function(x) data.frame(host = x[[2]]))
    host_mods <- dplyr::bind_rows(host_mods, .id = 'module')
    mods <- seq_along(mod_list)
  } else {
    if (!inherits(modules, 'data.frame') || !(c('name', 'module', 'type') %in% names(modules))) {
      stop('`modules` should be a `moduleWeb` or a data.frame with "name", "module" and "type" columns.')
    }
    host_mods <- modules %>%
      dplyr::filter(.data$type == "host") %>%
      dplyr::select(.data$name, .data$module) %>%
      dplyr::rename(host = .data$name)
    if (!is.null(module_order)) {
      mods <- module_order
    } else {
      mods <- sort(unique(host_mods$module))
    }
  }

  prepare_table <- function(m, s = NULL) {
    df <- as.data.frame(m)
    df$node <- rownames(df)
    df <- tidyr::pivot_longer(df, -.data$node, names_to = 'host')
    df <- df[df$value >= threshold, ]
    df$value <- NULL
    if (!is.null(s)) df$state <- s
    return(df)
  }

  # Get the ancestral states, reformat to plot on the parasite tree
  if (type == 'states') {
    l <- lapply(state, function(s) {
      prepare_table(at_nodes[['post_states']][, , as.character(s)], s)
    })
    node_df <- data.table::rbindlist(l)
  } else {
    node_df <- prepare_table(at_nodes[['post_repertoires']][, , repertoire])
  }

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
  suppressMessages( #suppress "scale for `y` is already present message
    p <- ggtree::ggtree(tree, layout = layout, ladderize = ladderize) +
      ggplot2::scale_x_continuous(name = NULL, labels = abs, expand = ggplot2::expansion(c(0.05, 0))) +
      ggplot2::scale_y_continuous(expand = c(0, 0.5)) +
      color_scale
  )
  # Flip the time axis the right way around
  if (layout %in% c('rectangular', 'slanted')) {
    p <- ggtree::revts(p + ggtree::theme_tree2())
  }

  # Extract the node coordinates, so we can easily plot our own node information
  coords <- p$data[, c('x', 'y', 'label', 'isTip')]
  coords <- coords[order(coords$y), ]

  x_range <- diff(range(coords$x))
  y_range <- diff(range(coords$y))
  node_df2 <- dplyr::left_join(node_df, coords, by = c('node' = 'label')) %>%
    dplyr::arrange(.data$node, .data$module) %>%
    dplyr::group_by(.data$node)

  if (layout %in% c('rectangular', 'slanted')) {
    node_df2 <- dplyr::mutate(node_df2, x = offset_x(.data$x, width = x_range * dodge_width))
    if (is.null(point_shape)) point_shape <- 15
  } else {
    node_df2 <- dplyr::mutate(node_df2, y = offset_x(.data$y, width = y_range * dodge_width))
    if (is.null(point_shape)) point_shape <- 16
  }

  if (type == "states" & length(state) > 1) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes_(~x, ~y, color = ~module, alpha = ~factor(.data$state, levels = .env$state)),
      node_df2, shape = point_shape, size = point_size
    ) +
      ggplot2::scale_alpha_ordinal(
        limits = factor(1:2), range = state_alpha, name = 'Interaction type',
        labels = c('1' = 'Potential', '2' = 'Actual'),
        guide = guide
      )

    if(threshold < 0.5) warning("Threshold is smaller than 0.5, so both states might be displayed for the same interaction.")

    } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes_(~x, ~y, color = ~module),
      node_df2, shape = point_shape, size = point_size
    )
  }
  return(p)
}

#' Plot a network with modules as an adjacency matrix, with aligned phylogenies.
#'
#' @param net An adjacency matrix for a bipartite network. This should be the extant network.
#'   Parasites should be the rows, hosts should be columns. If all values are 0 or 1 an binary
#'   network is represented, otherwise a weighted network is assumed.
#' @param tree The phylogeny of the symbiont clade (e.g. parasites, herbivores). Object of class `phylo`.
#' @param host_tree The phylogeny belonging to the hosts. Object of class `phylo`.
#'
#' See the examples on how to change the color scale.
#'
#' @return An assembly of plots, of class `patchwork`.
#' @export
#'
#' @inheritParams plot_ancestral_states
#' @inheritParams plot_extant_matrix
#'
#' @examples
#' \dontrun{
#'   at_nodes <- posterior_at_nodes(history, tree, host_tree, 66 + 1:65)
#'   plot_matrix_phylo(extant_net, at_nodes, tree, host_tree)
#'   # manual_colors
#'   plot_matrix_phylo(extant_net, at_nodes, tree, host_tree, colors = rainbow(20))
#' }
plot_matrix_phylo <- function(
    net, at_nodes, tree, host_tree, type = "states", state = 2, repertoire = 'fundamental',
    modules = NULL, module_order = NULL,
    threshold = 0.9, point_size = 3, dodge_width = 0.025, colors = NULL, ladderize = FALSE
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
  if (is.null(at_nodes)){
    p <- ggtree::ggtree(tree, ladderize = ladderize)
    parasite_plot <- ggtree::revts(p + ggtree::theme_tree2())
  } else {
    parasite_plot <- plot_ancestral_states(
      tree, at_nodes, modules, module_order, type = type, state = state,
      repertoire = repertoire, threshold = threshold, point_size = point_size,
      dodge_width = dodge_width, legend = FALSE, colors = colors, ladderize = ladderize
    )
  }

  parasite_coords <- parasite_plot$data[parasite_plot$data$isTip, c('x', 'y', 'label', 'isTip')]
  parasite_coords <- parasite_coords[order(parasite_coords$y), ]

  # Make the host tree plot
  if (inherits(modules, 'moduleWeb')) {
    mods <- seq_along(bipartite::listModuleInformation(modules)[[2]])
  } else {
    mods <- sort(unique(modules$module))
  }
  if (!is.null(module_order)) {
    mods <- module_order
  }

  host_plot <- ggplot2::ggplot(host_tree) +
    ggtree::geom_tree() +
    ggplot2::coord_flip() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0.5)) +
    ggtree::theme_tree()
  host_coords <- host_plot$data[host_plot$data$isTip, c('y', 'label')]
  host_coords <- host_coords[order(host_coords$y), ]

  # Make the matrix
  module_plot <- plot_extant_matrix(net, modules, module_order, parasite_coords$label, host_coords$label)
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



#' Plot ancestral networks with module information
#'
#' @param summary_networks A list of incidence matrices (summary network) for each time slice in
#'   `ages`. Output of `get_summary_network()`
#' @param matched_modules A list of lists containing the module information for each node at each
#'   network. Output of `modules_across_ages()`.
#' @param tree The phylogeny of the symbiont clade (e.g. parasites, herbivores). Object of class
#'   `phylo`.
#' @param module_levels Order in which the modules should be organized. Affects which color each
#'   module will be assigned. If NULL, takes the order of appearance in `matched_modules$module`.
#' @param palette Color palette used to plot module information.
#' @param node_size Size of the nodes in every network. If NULL, the default size is plotted.
#'
#' @return A list of plots of class `patchwork`. Each element contains the tree and network at a
#'   given time slice.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   data(tree)
#'   data(host_tree)
#'   data(history)
#'
#'   ages <- c(60, 50, 40, 0)
#'   at_ages <- posterior_at_ages(history, ages, tree, host_tree)
#'   summary_networks <- get_summary_networks(at_ages, threshold = 0.5, weighted = TRUE)
#'   all_mod <- modules_across_ages(summary_networks, tree)
#'
#'   plot <- plot_ancestral_networks(summary_networks, all_mod, tree)
#'   patchwork::wrap_plots(plot, guides = "collect")
#' }
plot_ancestral_networks <- function(
    summary_networks, matched_modules, tree, module_levels = NULL, palette = NULL, node_size = NULL
){

  # Input checking
  if (!is.list(summary_networks) || !all(vapply(summary_networks, inherits, TRUE, 'matrix'))) {
    stop('`summary_networks` should be a list of matrices, usually generated by `get_summary_networks`.')
  }
  if (!is.null(matched_modules) && (
    !inherits(matched_modules, 'list') && !inherits(matched_modules, 'data.frame')
  )) {
    stop('`matched_modules` should be of class `list` or `data.frame`.')
  }
  if (inherits(matched_modules, 'list')) {
    # extract the data.frame from match_modules that we need
    matched_modules <- matched_modules[['matched_modules']][['nodes_and_modules_per_age']]
  } else {
    if (inherits(matched_modules, 'data.frame') && !all(
      c('age','name', 'module') %in% names(matched_modules)
    )) {
      stop('`matched_modules` should be a `list` or a data.frame with "age", "name" and "module"columns.')
    }
  }
  if (!inherits(tree, 'phylo')) stop('`tree` should be of class `phylo`.')

  # ages has to be in decreasing order - '0' is the last element
  ages <- names(summary_networks) %>% as.numeric() %>% sort(decreasing = TRUE)

  # check that summary_networks is also in decreasing order
  # stop('`summary_networks` has to be in decreasing order of ages. Are you sure you used `get_summary_network`?'

  # check that all summary networks contain at least 3 rows (symbionts)
  # stop('all matrices in `summary_networks` have to contains at least 3 rows to produce a network')

  # auto-detect the network type
  net_vals <- sort(unique(unlist(summary_networks)))
  if (identical(net_vals, c(0, 1))) {
    weighted <- FALSE
    two_state <- FALSE
    weight_range <- c(0, 1)
  } else {
    if (identical(net_vals, c(0, 1, 2))) {
      weighted <- FALSE
      two_state <- TRUE
      weight_range <- c(0, 1)
    } else {
      weighted <- TRUE
      two_state <- FALSE
      weight_range <- range(net_vals[-1])
    }
  }
  # make list of tidy graphs
  list_tgraphs <- list()
  for(n in seq_along(summary_networks)){

    wnet <- summary_networks[[n]]

    graph <- tidygraph::as_tbl_graph(t(wnet), directed = F) %>%
      dplyr::left_join(dplyr::filter(dplyr::select(matched_modules, .data$age, .data$name, .data$module),
                                     .data$age == ages[n]), by = "name") %>%
      dplyr::select(.data$type, .data$name, .data$module)

    list_tgraphs[[n]] <- graph
  }

  # and root.time
  tree$root.time <- max(dispRity::tree.age(tree)$ages)

  # Slice the tree at ages and create data frame with module info
  list_subtrees <- list()
  list_tip_data <- list()

  # model "acctran" always uses the value from the ancestral node
  for(i in 1:(length(ages)-1)){
    subtree <- dispRity::slice.tree(tree, age = ages[[i]], "acctran")
    list_subtrees[[i]] <- subtree

    graph <- list_tgraphs[[i]]
    mod_from_graph <- tibble::tibble(
      module = tidygraph::activate(graph, what = nodes) %>%
        dplyr::filter(.data$type == TRUE) %>%
        dplyr::pull(.data$module),
      label = tidygraph::activate(graph, what = nodes) %>%
        dplyr::filter(.data$type == TRUE) %>%
        dplyr::pull(.data$name)
    )
    # extra step just to check that tip labels and graph node names match
    tip_data <- tibble::tibble(label = subtree$tip.label) %>%
      dplyr::inner_join(mod_from_graph, by = "label")
    list_tip_data[[i]] <- tip_data
  }

  # add tree and module info for present time
  list_subtrees[[length(ages)]] <- tree
  list_tip_data[[length(ages)]] <- tibble::tibble(label = tree$tip.label) %>%
    dplyr::inner_join(dplyr::filter(matched_modules, .data$age == 0),
                      by = c("label" = "name")) %>%
    dplyr::select(.data$label, .data$module)

  # plot
  if(is.null(module_levels)) module_levels <- sort(unique(matched_modules$module))

  plot_list <- list()

  for(t in 1:length(ages)){
    plot_age <- plot_network_at_age(
      list_subtrees[[t]], list_tip_data[[t]], list_tgraphs[[t]],
      module_levels, palette, tree, age = ages[t], weighted = weighted, two_state = two_state,
      weight_range = weight_range, node_size = node_size
    )
    plot_list[[t]] <- plot_age
  }

  return(plot_list)

}

nodes <- NULL


#' Plot one ancestral network with module information at a given time
#'
#' @param subtree a `phylo` object of the original tree sliced at a given time in the past.
#' @param tip_data a `data.frame` containing the module information for each tip in the subtree.
#' @param tgraph a `tbl_graph` containing the nodes and edges of the ancestral network and the
#'   module information for each node.
#' @param module_levels Order in which the modules should be organized. Affects which color each
#'   module will be assigned.
#' @param palette Color palette used to plot module information.
#' @param tree The phylogeny of the symbiont clade (e.g. parasites, herbivores). Object of class
#'   `phylo`.
#' @param age Age of the ancestral network to be plotted as the tittle.
#' @param weighted Whether the network should have weighted edges.
#' @param weight_range The range of weights conscidered for the width of the edges.
#' @param two_state Whether the width of the edges should reflect the state (instead of the
#'   posterior probability).
#'
#' @return An assembly of plots, of class `patchwork`.
#' @export
#' @importFrom ggtree %<+%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' plot_network_at_age(subtree, tip_data, tgraph, module_levels, palette, tree, age)
#' }
plot_network_at_age <- function(
    subtree, tip_data, tgraph, module_levels, palette = NULL, tree, age, weighted = TRUE,
    weight_range = c(0, 1), two_state = FALSE, node_size = NULL
) {

  if(is.null(palette)) palette <- scales::hue_pal()(length(module_levels))

  ggt <- ggtree::ggtree(subtree, ladderize = T, root.position = -tree$root.time) %<+% tip_data +
    ggtree::geom_tippoint(ggplot2::aes(color = factor(.data$module, levels = module_levels))) +
    ggtree::geom_rootedge(rootedge = 1) +
    ggplot2::scale_color_manual(values = palette, na.value = "grey70", drop = F) +
    ggtree::theme_tree2() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlim(c(-tree$root.time,0)) +
    ggplot2::labs(title = paste0(age, " Ma"))

  if (weighted) {
    edge_scale <- ggraph::scale_edge_width_continuous(limits = weight_range, range = c(0.2, 1))
    edge_aes <- ggplot2::aes(width = .data$weight)
  } else {
    if (two_state) {
      edge_scale <- ggraph::scale_edge_width_manual(
        name = 'State', guide = ggplot2::guide_legend(), values = c('1' = 0.2, '2' = 1),
        limits = c('1', '2')
      )
      edge_aes <- ggplot2::aes(width = as.character(.data$weight))
    } else {
      edge_scale <- NULL
      edge_aes <- NULL
    }
  }

  if(is.null(node_size)){
    geom_node_point <- ggraph::geom_node_point(
      ggplot2::aes(shape = .data$type, color = factor(.data$module, levels = module_levels))
    )
  } else {
    geom_node_point <- ggraph::geom_node_point(
      ggplot2::aes(shape = .data$type, color = factor(.data$module, levels = module_levels)),
      size = node_size
    )
  }

  ggn <- ggraph::ggraph(tgraph, layout = 'stress') +
    ggraph::geom_edge_link(edge_aes, color = "grey50") +
    geom_node_point +
    ggplot2::scale_shape_manual(
      values = c("square", "circle"), labels = c("Host", "Symbiont"), name = NULL
    ) +
    ggplot2::scale_color_manual(
      values = palette, na.value = "grey70", drop = F, name = "Module"#, limits = c(module_levels, NA)
    ) +
    edge_scale +
    ggplot2::theme_void()

  plot_age <- ggt + ggn + patchwork::plot_layout(widths = c(2,3))

  return(plot_age)

}
