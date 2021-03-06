# modified compute_modules() from bipartite
setClass("moduleWeb", representation(originalWeb="matrix", moduleWeb="matrix", orderA="vector", orderB="vector", modules="matrix", likelihood="numeric"));

#' Modified version of bipartite's computer_module function
#'
#' @param web Incidence matrix
#' @param method Becket
#' @param steps steps
#' @param tolerance tolerance
#' @param forceLPA FALSE
#'
#' @return Modularity Q
#' @export
#'
#' @examples
mycomputeModules = function(web, method="Beckett", steps=1000000, tolerance=1e-10, forceLPA=FALSE) {

  # check if, for binary data, any species is present everywhere ("empty" takes care of the "nowhere"):
  # if (length(table(unlist(web))) == 2  & ( any(colSums(web) == nrow(web)) | any(rowSums(web) == ncol(web)))){
  #   warning("Your binary data set contains one (or more) species present everywhere. These will be ignored, as they contain no information for the modularity algorithm.")
  #   nonWeb <- (!web) * 1
  #   web <- (! empty(nonWeb)) * 1
  # }

  # always uses DIRTLPA unless forced to use the faster LPA, which gets stuck in local optima more often (see Beckett 2016)
  web <- as.matrix(bipartite::empty(web))
  if (any(attr(web, "empty")) > 0) warning("Some empty columns or rows were deleted.")

  mod <- if (forceLPA) bipartite::LPA_wb_plus(web) else  bipartite::DIRT_LPA_wb_plus(web)
  # convert into moduleWeb-object:
  result <- bipartite::convert2moduleWeb(web,  mod)

  return(result)

}



# This function is ALSO in drawModules.R!! Auxiliary function checking whether
# the passed object is an object of class "moduleWeb" and contains correctly
# formatted information
isCorrectModuleWebObject = function(moduleWebObject) {

  if (!is(moduleWebObject, "moduleWeb")) {
    warning("Object of wrong class.");
    FALSE;
  }
  else if(dim(slot(moduleWebObject, "originalWeb")) == 0 ||  dim(slot(moduleWebObject, "moduleWeb")) != dim(slot(moduleWebObject, "originalWeb")) || dim(slot(moduleWebObject, "modules")) == 0) {
    warning("Object corrupt.");
    FALSE;
  }
  else if(min(slot(moduleWebObject, "originalWeb")) < 0 || min(slot(moduleWebObject, "moduleWeb")) < 0) {
    warning("entries of matrix have to be greater than or equal to 0.");
    FALSE;
  }
  else {
    TRUE;
  }
}
######----- ########

listModuleInformation = function(moduleWebObject) {

  if(isCorrectModuleWebObject(moduleWebObject)) {

    result	= list();

    web	= slot(moduleWebObject, "originalWeb");
    modules	= slot(moduleWebObject, "modules");

    n_a	= nrow(web);
    n_b	= ncol(web);

    offset	= 2;

    for(depth in unique(modules[,1])) {
      result[[depth+1]] = list();

      counter = 1;

      for(i in 1:nrow(modules)) {
        if(modules[i,1] == depth) {
          result[[depth+1]][[counter]]		= list();
          result[[depth+1]][[counter]][[1]]	= rownames(web)[modules[i,(offset+1):(n_a+offset)][modules[i,(offset+1):(n_a+offset)] > 0]];
          result[[depth+1]][[counter]][[2]]	= colnames(web)[(modules[i,(n_a+offset+1):(n_a+n_b+offset)][modules[i,(n_a+offset+1):(n_a+n_b+offset)] > 0])-n_a];

          counter = counter + 1;
        }
      }
    }

    result;
  }
}


printoutModuleInformation = function(moduleWebObject) {

  if(isCorrectModuleWebObject(moduleWebObject)) {

    modules	= slot(moduleWebObject, "modules");

    listOfModules = listModuleInformation(moduleWebObject);

    linebreak = "\n";

    a = linebreak;

    for(depth in unique(modules[,1])) {
      for(i in 1:length(listOfModules[[depth+1]])) {
        a = paste(a, "Depth: ", depth, linebreak, linebreak, "Nr of module: ", i, linebreak, linebreak, "Rownames: ", linebreak, sep="");
        for(j in 1:length(listOfModules[[depth+1]][[i]][[1]])) {
          a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[1]][j], sep=""), sep=linebreak);
        }
        a = paste(a, linebreak, linebreak, "Colnames: ", linebreak, sep="");
        for(j in 1:length(listOfModules[[depth+1]][[i]][[2]])) {
          a = paste(a, paste("\t", listOfModules[[depth+1]][[i]][[2]][j], sep=""), sep=linebreak);
        }
        a = paste(a, linebreak, linebreak, "__________________________", linebreak, linebreak, sep="");
      }
      a = paste(a, "__________________________", linebreak, linebreak, sep="");
    }

    cat(a);
  }
}

