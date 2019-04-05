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
    stop("This functionality is still under construction")
    # MIN-MAX NORMALIZATION
    hibound = c(22.2327575745000,118567.256530597,6.13594978300000,3.04444398050000,6.71199174400000,4.36844735975000,2.76309309050000)
    lobound = c(-11.5456138205000,-118558.932877979,-2.26230229700000,-0.740212034499999,-0.959183515999999,-1.03595686775000,-1.10192697450000)
    minX = c(1.16693895800000,-3.39917000000000e+16,0.0405368660000000,0.339137880000000,0,0,0)
    lambdaX = c(-0.223125000000000,0.528437500000001,1.30043750000000,-0.844812500000001,-0.794250000000001,1.91625000000000,0.830687500000001)
    muX = c(1.40873020972154,1001.03537420081,2.19669139716315,0.469945194965387,0.839479593949455,3.03822639738329,0.742915406888997)
    sigmaX = c(0.314519224026061,181.043417565308,0.858373279933531,0.112661854904320,0.0606316425681543,0.891291984344278,0.246925187328746)

  }else{
    # MIN-MAX & MEDIAN-IQR NORMALIZATION
    hibound = c(23.7389038738447,2.45595937752684,803.398085270330,2.43144861669896,2.88717834515708,15.8116522654644,0.340141846001221)
    lobound = c(-14.1225967541981,-0.709874366341159,-771.931418603663,-1.67225739529610,-0.785763036773620,-11.5751116064500,-0.277641846001221);
    minX = c(1.62357386179707,0,0.769230769230769,0,0.237625121348345,1.00567121540896,0);
    lambdaX =c(-0.205750000000000,0.501187500000001,-0.128250000000000,-1.75981250000000,-1.51662500000000,-0.665562500000001,-13.1663125000000);
    muX = c(1.22159342595184,0.715972404772611,2.43996997324290,0.236704355099472,0.404112336345564,0.623394589468880,0.0287167367721530);
    sigmaX = c(0.451052410971313,0.227282405566549,1.31057342664720,0.131392003960080,0.0703408312315028,0.360541179317276,0.0234034120625995);

  }
  for(i in 1:dim(x)[2]){
    recs <- which(x[ ,i] > hibound[i])
    x[recs, i] <- hibound[i]
    recs <- which(x[ ,i] < lobound[i])
    x[recs, i] <- lobound[i]

    x[ ,i] <- x[ ,i] - minX[i] + 1
    x[ ,i] <- (x[ ,i]^lambdaX[i] - 1)/lambdaX[i]
    x[ ,i] <- (x[ ,i] - muX[i])/sigmaX[i]

  }
  return(x)
}
