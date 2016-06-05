single_predict <- function(rit,instance,radius,absence,compute_counts){
  class_names <- names(rit[["Class_priors"]])
  nb_class <- length(class_names)
  
  probs <- vector(mode="numeric",length=nb_class) + log10(rit[["Class_priors"]])
  names(probs) <- class_names
  
  model <- rit[["Model"]]
  map <- rit[["Map"]]
  
  count_has <- 0
  for(elem in model){
    weight <- compute_weight(instance,elem[["interaction"]],map,radius)
    probs <- probs + weight*log10(elem[["sm"]])
    if(absence)
      probs <- probs + (1-weight)*log10(elem[["sm_neg"]])
    if(weight == 1)
      count_has <- count_has + 1
  }
  
  # Compute argmax
  response <- names(probs)[which.max(probs)]
  if(compute_counts){
    return(list(response=response,count_has=count_has/length(model)))
  }
  else{
    return(response)
  }
}

naive_predict <- function(rit,testset,radius=0,absence=FALSE,compute_counts=FALSE){
  if(radius < 0)
    stop("Radius should be positive.")
  
  if(compute_counts)
    response <- vector(mode="list",length=nrow(testset))
  else
    response <- vector(mode="character",length=nrow(testset))
  
  for(i in 1:nrow(testset)){
    response[[i]] <- single_predict(rit,testset[i,],radius,absence,compute_counts)
  }
  
  if(compute_counts){
    r <- sapply(response,function(i){
      i[["response"]]
    })
    c <- sapply(response,function(i){
      i[["count_has"]]
    })
    return(list(response=r,count_has=sum(c)))
  }
  else{
    return(response)
  }
}