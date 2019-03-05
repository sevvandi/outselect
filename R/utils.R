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


renormalize <- function(x, dd=2){

  if(dd==1){
    stop("This will be implemented later!")
  }else{
    hibound = c(23.7389038738447,2.45595937752684,803.398085270330,2.43144861669896,2.88717834515708,15.8116522654644,0.340141846001221)
    lobound = c(-14.1225967541981,-0.709874366341159,-771.931418603663,-1.67225739529610,-0.785763036773620,-11.5751116064500,-0.277641846001221);
    minX = c(1.62357386179707,0,0.769230769230769,0,0.237625121348345,1.00567121540896,0);
    lambdaX =c(-0.205750000000000,0.501187500000001,-0.128250000000000,-1.75981250000000,-1.51662500000000,-0.665562500000001,-13.1663125000000);
    muX = c(1.22159342595184,0.715972404772611,2.43996997324290,0.236704355099472,0.404112336345564,0.623394589468880,0.0287167367721530);
    sigmaX = c(0.451052410971313,0.227282405566549,1.31057342664720,0.131392003960080,0.0703408312315028,0.360541179317276,0.0234034120625995);


    for(i in 1:dim(x)[2]){
      recs <- which(x[ ,i] > hibound[i])
      x[recs, i] <- hibound[i]
      recs <- which(x[ ,i] < lobound[i])
      x[recs, i] <- lobound[i]

      x[ ,i] <- x[ ,i] - minX[i] + 1
      x[ ,i] <- (x[ ,i]^lambdaX[i] - 1)/lambdaX[i]
      x[ ,i] <- (x[ ,i] - muX[i])/sigmaX[i]

    }
  }

  return(x)
}
