#'Performs \code{binom.test} to see if Min-Max is actually better than Median-IQR method.
#'
#'@param rocpr If \code{rocpr=1} then area under ROC curve is used as the performance measure.  If \code{rocpr=2} area under PR curve is used.
#'@param prop Claimed proportion that Min-Max is better. Used for the null hypothesis.
#'
#'@return A list containing the following:
#' \describe{
#'  \item{percentages}{The proportion of datasets that gave better performance for Min-Max.}
#'  \item{pvalues}{The p-value of the  \code{binom.test}.}
#'  \item{confintervals}{The confidence intervals of the  \code{binom.test}.}
#' }
#'
#' @examples
#' out <- IsMinMaxBetter(1)
#' out$confintervals
#' @export

IsMinMaxBetter <- function(rocpr=1, prop=0.5){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }
  e <- new.env()
  if(rocpr==1){
    # ROC values
    data("perf_vals_roc_all", envir=e)
    perfs <- perf_vals_roc_all[ ,2*(1:(dim(perf_vals_roc_all)[2]/2))]

  }else{
    # PR values
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all[ ,2*(1:(dim(perf_vals_pr_all)[2]/2))]
  }
  perfs_mm <- perfs[ , 2*(1:(dim(perfs)[2]/2))]
  perfs_iq <- perfs[ , 2*(1:(dim(perfs)[2]/2))-1]

  perfs_diff <- perfs_mm - perfs_iq
  percentages <- apply(perfs_diff,2,function(x) sum(x>0)/length(x) )
  p_values <- rep(0, dim(perfs_diff)[2])
  conf_ints <- matrix(0,nrow=dim(perfs_diff)[2], ncol=2 )
  for(i in 1:dim(perfs_diff)[2]){
    binom_test <- binom.test(sum(perfs_diff[ ,i]>0),dim(perfs_diff)[1],prop, alternative="less", conf.level = 0.99)
    p_values[i] <- binom_test$p.value
    conf_ints[i, ] <- binom_test$conf.int[1:2]
  }
  out <- list()
  out$percentages <- percentages
  out$pvalues <- p_values
  out$confintervals <- conf_ints
  return(out)
}

#' Computes percentages of datasets which are sensitive to normalization.
#' @inheritParams IsMinMaxBetter
#' @param xi Sensitivity to normalization parameter. 	For a given dataset, we say that an outlier detection method is \code{xi}-sensitive to normalization if the difference between the maximum performance and the minimum performance across all normalization schemes for that outlier detection method is greater than \code{xi}.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{percentages}{The proportion of datasets that gave better performance for Min-Max.}
#'  \item{pvalues}{The p-value of the  \code{binom.test}.}
#'  \item{confintervals}{The confidence intervals of the  \code{binom.test}.}
#'  \item{methods}{The outlier detection methods.}
#' }
#'
#'@examples
#'out <- SensitivityToNorm(1,0.05)
#'out$confintervals
#'
#'@importFrom stats binom.test
#'
#'@export

SensitivityToNorm <- function(rocpr, xi){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }
  if((xi>=1)|(xi<=0)){
    stop("Invalid xi. xi need to be between 0 and 1.")
  }
  e <- new.env()
  if(rocpr==1){
    # ROC values
    data("perf_vals_roc_all", envir=e)
    perfs <- perf_vals_roc_all
  }else{
    # PR values
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all
  }
  st_col <- 1
  en_col <- 4
  num_methods <- dim(perfs)[2]/4

  percentages <- rep(0,num_methods)
  p_values <- rep(0, num_methods)
  conf_ints <- matrix(0,nrow=num_methods, ncol=2 )
  methods <- c()

  for(i in 1:num_methods){
    perf_method <- perfs[ ,st_col:en_col]
    rangediff <- apply(perf_method, 1, function(x) diff(range(x, na.rm=TRUE)) )
    percentages[i] <- sum(rangediff > xi )/length(rangediff)

    binom_test <- binom.test(sum(rangediff > xi ),length(rangediff), 0.5, alternative="two.sided", conf.level = 0.99)
    p_values[i] <- binom_test$p.value
    conf_ints[i, ] <- binom_test$conf.int[1:2]

    pos <- regexpr("_", colnames(perfs)[st_col])[1] -1
    method <- substring(colnames(perfs)[st_col], 1, pos)
    methods <- c(methods, method)

    st_col <- en_col + 1
    en_col <- en_col + 4
  }
  methods[which(methods=="FAST")] <- "FAST_ABOD"
  out <- list()
  out$percentages <- percentages
  out$pvalues <- p_values
  out$confintervals <- conf_ints
  out$methods <- methods
  return(out)
}
