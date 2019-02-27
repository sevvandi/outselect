#'Performs \code{binom.test} to see if Min-Max is actually better than Median-IQR method.
#'
#'@param rocpr If \code{rocpr=1} then area under ROC curve is used as the performance measure.  If \code{rocpr=2} area under PR curve is used.
#'@param prop Claimed proportion that Min-Max is better. Used for the null hypothesis.
#'@param tt If \code{tt=1}, then the two sided binomial test is performed. If \code{tt=2} then the alternative is less than.
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

IsMinMaxBetter <- function(rocpr=1, prop=0.5, tt=1){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }
  if((tt!=1)&(tt!=2)){
    stop("Invalid tt. tt should equal 1 or 2.")
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
    binom_test <- binom.test(sum(perfs_diff[ ,i]>0),dim(perfs_diff)[1],prop, alternative=c("two.sided", "less")[tt], conf.level = 0.99)
    p_values[i] <- binom_test$p.value
    conf_ints[i, ] <- binom_test$conf.int[1:2]
  }
  out <- list()
  out$percentages <- percentages
  out$pvalues <- p_values
  out$confintervals <- conf_ints
  return(out)
}


#' Computes sensitive to normalization statistics.
#' @inheritParams IsMinMaxBetter

#' @return A list containing the following:
#' \describe{
#'  \item{friedman}{The output of the Friedman test.}
#'  \item{nemenyi}{The output of the Nemenyi test \code{tsutils::nemenyi}.}
#' }
#'
#'
#'
#'
#'@examples
#'out <- SensitivityToNorm(1)
#'out
#'
#'@importFrom stats binom.test
#'
#'@export

SensitivityToNorm <- function(rocpr=1){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }

  e <- new.env()
  if(rocpr==1){
    # ROC values
    data("perf_vals_roc_all", envir=e)
    perfs <- perf_vals_roc_all
    data("filenames_roc", envir=e)
    filenames <- filenames_roc
  }else{
    # PR values
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all
    data("filenames_pr", envir=e)
    filenames <- filenames_pr
  }
  st_col <- 1
  en_col <- 4
  num_methods <- dim(perfs)[2]/4

  # percentages <- rep(0,num_methods)
  # conf_ints <- matrix(0,nrow=num_methods, ncol=2 )
  # means <- rep(0, num_methods)
  # sds <- rep(0, num_methods)
  # rates <- rep(0,num_methods)
  rangediff <- matrix(0, nrow=dim(perfs)[1], ncol=num_methods)
  methods <- c()


  # Max - Min for each outlier method across the normalization methods
  for(i in 1:num_methods){
    perf_method <- perfs[ ,st_col:en_col]
    rangediff[,i] <- apply(perf_method, 1, function(x) diff(range(x, na.rm=TRUE)) )

    pos <- regexpr("_", colnames(perfs)[st_col])[1] -1
    method <- substring(colnames(perfs)[st_col], 1, pos)
    methods <- c(methods, method)

    st_col <- en_col + 1
    en_col <- en_col + 4
  }

  # ---- Find the sources for filenames
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
  uniq_f_s <- unique(file_source)
  methods[which(methods=="FAST")] <- "FAST_ABOD"
  colnames(rangediff) <- methods

  # medians <- apply(rangediff, 2, median)
  # ordering <- order(medians)
  # rangediff <- rangediff[, ordering]
  # methods <- methods[ordering]

  df <- cbind.data.frame(file_source, rangediff)
  # --- Make a big data frame with source, outlier method and sensitivity to normalization
  for(j in 1:num_methods){
    temp <- reshape::melt(df[ ,c(1,(j+1))] )
    colnames(temp) <- c("source", "method", "value")
    if(j==1){
      dat <- temp
    }else{
      dat <- rbind.data.frame(dat, temp)
    }
  }

  df2 <- aggregate(dat$value, by=list(m=dat$method, s=dat$source), FUN=median)

  friedman_test <- stats::friedman.test(y=df2$x, groups=df2$m, blocks=df2$s)
  # nemenyi <- PMCMR::posthoc.friedman.nemenyi.test(y=df2$x, groups=df2$m, blocks=df2$s)

  df3 <- aggregate(df[,-1], by=list(file_source), FUN=median)
  nemenyi <- tsutils::nemenyi(as.matrix(df3[ ,-1]), conf.level=0.95, sort=TRUE, plottype="vline")

  out <- list()
  out$friedman <-friedman_test
  out$nemenyi <- nemenyi
  return(out)
}




SensitivityToNormMixedMod <- function(rocpr=1){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }

  e <- new.env()
  if(rocpr==1){
    # ROC values
    data("perf_vals_roc_all", envir=e)
    perfs <- perf_vals_roc_all
    data("filenames_roc", envir=e)
    filenames <- filenames_roc
  }else{
    # PR values
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all
    data("filenames_pr", envir=e)
    filenames <- filenames_pr
  }

  # ---- Find the sources for filenames
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
  uniq_f_s <- unique(file_source)

  st_col <- 1
  en_col <- 4
  num_methods <- dim(perfs)[2]/4
  rangediff <- matrix(0, nrow=dim(perfs)[1], ncol=num_methods)
  methods <- c()
  for(i in 1:num_methods){
    perf_method <- perfs[ ,st_col:en_col]
    rangediff[,i] <- apply(perf_method, 1, function(x) diff(range(x, na.rm=TRUE)) )

    pos <- regexpr("_", colnames(perfs)[st_col])[1] -1
    method <- substring(colnames(perfs)[st_col], 1, pos)
    methods <- c(methods, method)
    st_col <- st_col + 4
    en_col <- en_col + 4
  }
  methods[which(methods=="FAST")] <- "FAST_ABOD"
  colnames(rangediff) <- methods
  df <- cbind.data.frame(file_source, rangediff)

  # --- Make a big data frame with source, outlier method and sensitivity to normalization
  for(j in 1:num_methods){
    temp <- reshape::melt(df[ ,c(1,(j+1))] )
    colnames(temp) <- c("source", "method", "value")
    if(j==1){
      dat <- temp
    }else{
      dat <- rbind.data.frame(dat, temp)
    }
  }

  fit <- lme4::lmer(value ~  method + (1 | source), data=dat)
  obj <- multcomp::glht(fit, linfct = multcomp::mcp(method = "Tukey"))

  out <- list()
  out$fit <- fit
  out$glht <- obj
  return(out)

}
