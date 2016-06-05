single_predict <- function(rit,instance){
  class_names <- names(rit[["Class_priors"]])
  nb_class <- length(class_names)
  
  probs <- vector(mode="numeric",length=nb_class) + log10(rit[["Class_priors"]])
  names(probs) <- class_names
  
  model <- rit[["Model"]]
  map <- rit[["Map"]]
  
  for(elem in model){
    weight <- compute_weight(instance,elem[["interaction"]],map)
    probs <- probs + weight*log10(elem[["sm"]])
  }
  
  # Compute argmax
  response <- names(probs)[which.max(probs)]
  response
}

#' @title Classification Rule for Random Intersection Trees with discretization.
#' @description Applies a basic \code{argmax} rule in order to classify new instances.
#'
#' @return A response vector for the \code{testset} instances
#'
#' @param rit A model produced by \code{naive_RIT}
#' @param testset A dataframe containing the instance to classify
#' 
#' @references Ballarini Robin. Random intersection trees for genomic data analysis. Master's thesis, Université Catholique de Louvain, 2016.
#' @export
#'
naive_predict <- function(rit,testset){
  response <- vector(mode="character",length=nrow(testset))
  
  for(i in 1:nrow(testset)){
    response[[i]] <- single_predict(rit,testset[i,])
  }
  
  response
}