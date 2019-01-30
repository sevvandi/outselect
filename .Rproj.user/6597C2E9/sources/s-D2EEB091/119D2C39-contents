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
