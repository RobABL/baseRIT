# Equal width unsupervised discretization
equal_width <- function(data,n=nrow(data)^(1/3)){
  cutp <- lapply(data,function(attr){
    seq(min(attr),max(attr),length.out=n+1)
  })
  
  disc <- sapply(1:ncol(data),function(attr_idx){
    findInterval(data[[attr_idx]],cutp[[attr_idx]],all.inside=T)
  })
  disc <- as.data.frame(disc)
  names(disc) <- names(data)
  
  list(cutp=cutp,Disc.data=disc)
}

naive_preproc <- function(data,classes,algo){
  # Empty dataframe with same number of rows
  new_data <- data.frame("rows"=1:nrow(data), row.names="rows")
  
  # Mapping between original and discretized attributes, used for predicting new instances
  map <- vector("list",length=ncol(data))
  
  # Count discretized attributes added
  count <- 1
  
  # Turn ordered factors into numeric vectors
  data <- lapply(data,function(col){
    if(is.ordered(col)){
      as.numeric(col)
    }
    else{
      col
    }
  })
  data <- as.data.frame(data)
  
  # Indexes
  cat <- vapply(data,is.factor,FALSE,USE.NAMES=FALSE)
  cont_idx <- which(!cat)
  cat_idx <- which(cat)
  
  # Continous attributes
  if(length(cont_idx) > 0){
    # Restrict df to continuous attributes + disc
    cont_data <- data[cont_idx]
    cont_vars <- sapply(cont_data,var,USE.NAMES=FALSE)
    
    if(algo == "caim"){
      cont_data[["Class"]] <- classes
      disc_out <- disc.Topdown(cont_data,method=1) # CAIM algorithm
    }
    else if(algo == "ameva"){
      cont_data[["Class"]] <- classes
      disc_out <- disc.Topdown(cont_data,method=3) # AMEVA algorithm
    }
    else if(algo == "chim"){ # chi-merge with level of significance 0.05
      cont_data[["Class"]] <- classes
      disc_out <- chiM(cont_data,alpha=0.05)
    }
    else if(algo == "chim1"){ # chi-merge with level of significance 0.1
      cont_data[["Class"]] <- classes
      disc_out <- chiM(cont_data,alpha=0.1)
    }
    else if(algo == "eqwidth"){ # equal-width binning (unsupervised)
      disc_out <- equal_width(cont_data)
    }
  
    # Add disc continuous attributes
    for(i in 1:length(cont_idx)){
      cuts <- disc_out$cutp[[i]]
      
      if(algo == "chim" || algo == "chim1"){
        cuts <- c(-Inf,cuts,Inf)
      }
      else{
        cuts[1] <- -Inf
        cuts[length(cuts)] <- Inf
      }
      
      for(j in 1:length(cuts)){
        if(j < length(cuts)){
          new_data[[count]] <- as.numeric(disc_out$Disc.data[[i]] == j)
          map[[count]] <- list(attr=cont_idx[i], value=c(cuts[j],cuts[j+1]),var=cont_vars[i])
          count <- count + 1
        }
      }
    }
  }
  
  # Cat attributes
  if(length(cat_idx) > 0){
    cat_data <- data[cat_idx]
    for(i in 1:length(cat_idx)){
      nb_values <- nlevels(cat_data[[i]])
      for(j in 1:nb_values){
        new_data[[count]] <- as.numeric(as.numeric(cat_data[[i]]) == j)
        map[[count]] <- list(attr=cat_idx[i], value=j)
        count <- count + 1
      }
    }
  }
    
  list(data=new_data,map=map)
}

naive_RIT <- function(data,classes,theta,disc="ameva",branch=3,depth=10L,split_nb=1,n_trees=100L,
                      min_inter_sz=2L,L=100L,es=TRUE){
  # check parameters
  if(!is.data.frame(data))
    stop("Data parameter must be a dataframe.")
  if(length(classes) != nrow(data))
    stop("Data and response vector must have the same number of rows.")
  if(nrow(data) == 0)
    stop("Data cannot be empty.")
  if(any(is.na(classes)))
    stop("Response vector cannot contain missing values.")
  if(!(disc %in% c("caim","ameva","chim","chim1","eqwidth")))
    stop("Unknown discretization algorithm.")
  if(split_nb < 1)
    stop("split_nb should be at least 1.")
  
  classes <- factor(classes)
  class_names <- levels(classes)
  nb_class <- nlevels(classes)
  
  if(missing(theta)){
    theta <- vector(mode="numeric",length=nb_class) - 1
    names(theta) <- class_names
  }
  else{
    if(length(theta) == 1){
      theta <- vector(mode="numeric",length=nb_class) + theta
      names(theta) <- class_names
    }
    else if(!(all(class_names %in% names(theta)))){
      stop("Named vector theta doesn't contain prevalence thresholds for some classes.")
    }
  }
  
  preproc <- naive_preproc(data,classes,disc)
  data <- preproc$data
  map <- preproc$map
  
  # Split dataset per class and order theta vector
  datas <- vector("list",length=nb_class)
  names(datas) <- class_names
  o_theta <- vector(mode="numeric",length=nb_class)
  names(o_theta) <- class_names
  for(cls in class_names){
    datas[[cls]] <- data[which(classes==cls),]
    o_theta[[cls]] <- theta[[cls]]
  }
  
  # Call main algo
  all_leaves <- cpp_naive_RIT(datas,o_theta,n_trees,depth,split_nb,branch,min_inter_sz,L,es)
  
  # Compute class prior probabilities
  priors <- data.table(classes)[,.N,keyby=classes]
  priors[,"N"] <- priors[,N] / length(classes)
  priors_vec <- vector(mode="numeric",length=nb_class)
  names(priors_vec) <- class_names
  for(i in 1:nrow(priors)){
    priors_vec[[priors[i,classes]]] <- priors[i,N]
  }
  
  # Smoothing
  cls_sums <- sapply(datas,nrow)
  sm_leaves <- lapply(all_leaves,function(inter){
    nb <- inter[["prevalence"]]*cls_sums
    nb_neg <- {1-inter[["prevalence"]]}*cls_sums
    list(interaction=inter[["interaction"]],prevalence=inter[["prevalence"]],sm={nb+1}/{nb+nb_neg+2},sm_neg={nb_neg+1}/{nb+nb_neg+2},depth=inter[["depth"]])
  })
  
  list(Model=sm_leaves,Map=map,Class_priors=priors_vec) 
}