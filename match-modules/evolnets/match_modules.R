#' Match modules across time slices
#'
#' @param module_info
#' @param ages
#' @param tree
#'
#' @return
#' @export
#'
#' @examples
match_modules <- function(module_info, ages, tree){

  # Make sure that ages are ordered present -> past and start at present
  ages <- sort(ages)
  if(ages[1] > 0) ages <- c(0,ages)

  # Prepare for adding children module information
  mod_df <- mutate(module_info,
                   child1_mod = 0,
                   child2_mod = 0,
                   module = case_when(age == 0 ~ original_module))

  # Get a tibble with tree information for finding child nodes
  tree_t <- as_tibble(tree)

  # Get subtree at each age
  tree$root.time <- max(tree.age(tree)$ages)
  list_subtrees <- list()
  for(i in 1:length(ages)){

    subtree <- slice.tree(tree, age = ages[[i]], "acctran")
    list_subtrees[[i]] <- subtree
  }

  # Find if tip is present in the next age and, if not, get children nodes
  for(t in 2:length(ages)){

    subtree <- list_subtrees[[t]]
    tips <- subtree$tip.label

    # look at next time step
    tips_next <- list_subtrees[[t-1]]$tip.label

    for(tip in 1:length(tips)){

      # if branch doesn't split until next time step
      if(tips[tip] %in% tips_next){

        # get module at the next time step
        mod_next <- module_info %>%
          filter(age == ages[t-1],
                 name == tips[tip]) %>%
          pull(original_module)

        # add that info to the data frame
        mod_df[which(mod_df$name == tips[tip] & mod_df$age == ages[t]),
               'child1_mod'] = mod_next

      } else{
        node_treeio <- tree_t %>%
          filter(label == tips[tip]) %>%
          pull(node)

        children <- tree_t %>%
          filter(parent == node_treeio) %>%
          pull(label)

        # find modules of both children
        mod_children <- module_info %>%
          filter(age == ages[t-1],
                 name %in% children) %>%
          pull(original_module)

        # get the size of the modules (the number of within-module interactions)
        # same thing for each child

        # get the geometric mean of within-module interactions
        # same thing for each child

        # add that info to the data frame
        mod_df[which(mod_df$name == tips[tip] & mod_df$age == ages[t]),
               'child1_mod'] <- mod_children[1]
        mod_df[which(mod_df$name == tips[tip] & mod_df$age == ages[t]),
               'child2_mod'] <- mod_children[2]

      }

    }

  }

  #names(list_subtrees) <- ages

  return(mod_df)

}



