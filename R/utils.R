EpsilonGood <- function(x, eps){
  ## x contains 12 columns for all datasets
  ## eg. eps = 0.05
  max_perf <- apply(x,1, function(x) max(x, na.rm = TRUE))
  output <- matrix(0,nrow=dim(x)[1], ncol=dim(x)[2])
  colnames(output) <- colnames(x)
  for(ii in 1:dim(x)[2]){
    perf.vals <- x[,ii]
    output[,ii] <- ifelse(perf.vals >=(max_perf-eps),1,0 )
  }
  return(output)
}


RemoveInfiniteValues <- function(X, cols){
  # X has the features
  # Cols are columns with infinite values
  iqr.vals <- apply(X[,cols],2, function(x)IQR(x, na.rm=TRUE))
  for(ii in 1:length(cols)){
    recs<- which(is.infinite(X[,cols[ii]]))
    X[recs,cols[ii]] <- 5*iqr.vals[ii]*sign(X[recs,cols[ii]])
  }
  return(X)
}

GetFileSources <- function(filenames){
  file_source <-c()
  for(ll in 1:length(filenames)){
    fname <- filenames[ll]
    regobj1 <- regexpr("_C", fname)
    regobj2 <- regexpr("_withoutdupl", fname)
    if(regobj1[1]<0){
      regobj <- regobj2
    }else if(regobj2[1]<0){
      regobj <- regobj1
    }else{
      regobj <- regobj1
    }
    end.ind <- regobj[1]-1
    file_source <- c(file_source, substring(fname, 1, end.ind))
  }
  return(file_source)
}
