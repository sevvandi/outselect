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
#'  \item{rangediff}{For each outlier method and each dataset, the maximum performance minus the minimum performance. }
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
  file_source <- GetFileSources(filenames)

  uniq_f_s <- unique(file_source)
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

  df2 <- stats::aggregate(dat$value, by=list(m=dat$method, s=dat$source), FUN=median)

  friedman_test <- stats::friedman.test(y=df2$x, groups=df2$m, blocks=df2$s)
  # nemenyi <- PMCMR::posthoc.friedman.nemenyi.test(y=df2$x, groups=df2$m, blocks=df2$s)

  df3 <- stats::aggregate(df[,-1], by=list(file_source), FUN=median)
  nemenyi <- tsutils::nemenyi(as.matrix(df3[ ,-1]), conf.level=0.95, sort=TRUE, plottype="vline", main="Nemenyi test average ranks")

  out <- list()
  out$friedman <-friedman_test
  out$nemenyi <- nemenyi
  out$rangediff <- rangediff
  return(out)
}


#' Computes sensitive to normalization for dataset sources
#' @inheritParams IsMinMaxBetter
#' @param xi The xi-sensitivity to normalization parameter xi. Defaults to \code{0.05}
#'
#'
#' @return A list containing the following:
#' \describe{
#'  \item{KruskalWallis}{The output of the KruskalWallis from \code{stats::kruskal.test}.}
#'  \item{max_val}{The dataset source, which is most sensitive to normalization.}
#'  \item{min_val}{The dataset source, which is least sensitive to normalization. }
#'  \item{xi_sensitive_num}{For how many outlier methods each dataset is sensitive to normalization.}
#' }
#'
#'
#'
#'
#'@examples
#'out <- SensitivityNormDatasets(1, 0.05)
#'out
#'
#'@export

SensitivityNormDatasets <- function(rocpr=1, xi=0.05){
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
  file_source <- GetFileSources(filenames)
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

  xi_sensitive_num <- apply(rangediff, 1, function(x)sum(x > xi))
  df <- cbind.data.frame(file_source, xi_sensitive_num)
  krusk <- stats::kruskal.test(xi_sensitive_num ~ file_source, data=df)

  df2 <- stats::aggregate(df[,-1], by=list(file_source), FUN=mean)
  range_min <- df2[which.min(df2[,2]), ]
  range_max <- df2[which.max(df2[,2]), ]
  colnames(df2) <- c("file_source", "Num_Out_Methods_Sensitive")

  out <- list()
  out$KruskalWallis <- krusk
  out$max_val <- range_max
  out$min_val <- range_min
  out$xi_sensitive_num <- xi_sensitive_num

  return(out)
}


#' Fits a mixed effects model to account for normalization methods and outlier methods accounting for the dataset variants
#' @inheritParams IsMinMaxBetter
#'
#' @return A list containing the following:
#' \describe{
#'  \item{fit1}{The first model without interactions between normalization and outlier method.}
#'  \item{fit2}{The second model with interactions between normalization and outlier method.}
#'  \item{aov}{The output of an ANOVA test. }
#' }
#'
#'
#'
#'
#'@examples
#'out <- SensitivityToNormMixedMod(1)
#'out
#'
#'@export


SensitivityToNormMixedMod <- function(rocpr=1, ds=FALSE, pp=FALSE){
  # if rocpr =1 then it is roc values, if 2 it is pr values
  # ds is a switch for down sampling. If ds == TRUE, then we down sample the observations such that each dataset source has a similar number of dataset variants generated.
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }
  print("This will take some time. . . . ")

  e <- new.env()
  if(rocpr==1){
    # ROC values
    data("perf_vals_roc_all", envir=e)
    data("filenames_roc", envir=e)
    df_normed <- cbind.data.frame(e$filenames_roc, e$perf_vals_roc_all)
    colnames(df_normed) <- c("filename", colnames(e$perf_vals_roc_all))

    data("not_normed_perfs_roc_14", envir=e)
    df_ori <- e$not_normed_perfs_roc_14

    df <- merge(df_ori, df_normed)


    filenames <- df$filename
    perfs <- df[ ,-1]
  }else{
    # PR values
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all
    data("filenames_pr", envir=e)
    filenames <- filenames_pr
  }

  # ---- Find the sources for filenames
  file_source <-GetFileSources(filenames)
  uniq_f_s <- unique(file_source)

  if(ds){
    # Downsampling is set to TRUE
    # Down sample such that each source dataset has approximately a fixed number of variants
    set.seed(1)
    pereach <- 10
    tab_fs <- table(file_source)
    rec_sel <- c()
    for(l in 1:length(uniq_f_s)){
     if(tab_fs[l]>pereach){
       rec_fs <- which(file_source %in% names(tab_fs)[l])
       temp <- sample(rec_fs, pereach)
       if(l==1){
         rec_sel <-  temp
       }else{
         rec_sel <- c(rec_sel, temp)
       }
     }
    }
    file_source <- file_source[rec_sel]
    uniq_f_s <- unique(file_source)
    filenames <- filenames[rec_sel]
    perfs <- perfs[rec_sel, ]
  }
  # ---- Norm method and outlier method
  norm.method <- c()
  outlier.method <- c()
  for(ii in 1:dim(perfs)[2]){
    cname2 <- colnames(perfs)[ii]
    regobj0 <- regexpr("_Original", cname2)
    regobj <- regexpr("_M", cname2)
    regobj1 <- ifelse(regobj0[1]>0, regobj0, regobj)
    out.m <- substring(cname2, 1,(regobj1[1]-1))
    norm.m <- substring(cname2, (regobj1[1]+1),nchar(cname2))
    norm.method <- c(norm.method,norm.m)
    outlier.method <- c(outlier.method,out.m )
  }

  # ---- Make a long data frame
  num_recs <- dim(perfs)[1]
  for(kk in 1:dim(perfs)[2]){
    cname2 <- colnames(perfs)[kk]
    regobj1 <- regexpr("_M", cname2)
    regobj0 <- regexpr("_Original", cname2)
    regobj <- regexpr("_M", cname2)
    regobj1 <- ifelse(regobj0[1]>0, regobj0, regobj)

    out.m <- substring(cname2, 1,(regobj1[1]-1))
    norm.m <- substring(cname2, (regobj1[1]+1),nchar(cname2))

    if(kk==1){
      dat.long <- cbind.data.frame(file_source,   rep(out.m, num_recs), rep(norm.m, num_recs), perfs[,kk] )
    }else{
      temp <-  cbind.data.frame(file_source,  rep(out.m,num_recs), rep(norm.m,num_recs), perfs[,kk] )
      dat.long <- rbind.data.frame(dat.long,temp)
    }

  }

  colnames(dat.long) <- c("Source","Out", "Norm",  "performance")
  fit.1 <- lme4::lmer(performance ~ Norm + Out + (1 | Source), data=dat.long)

  if(rocpr==1){
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Original"] <- "O"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Mean_SD"] <- "D"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Median_IQR"] <- "Q"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Median_MAD"] <- "M"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Min_Max"] <- "X"
    levels(dat.long$Norm) <- c("O", "D", "Q", "M", "X")
  }else{
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Mean_SD"] <- "D"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Median_IQR"] <- "Q"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Median_MAD"] <- "M"
    levels(dat.long$Norm)[levels(dat.long$Norm) == "Min_Max"] <- "X"
    levels(dat.long$Norm) <- c("D", "Q", "M", "X")
  }



  levels(dat.long$Out)[3] <- "F.ABOD"
  levels(dat.long$Out)[2] <- "Ens"

  fit.2 <- lme4::lmer(performance ~ Out*Norm + (1 | Source), data=dat.long)
  aov_obj <- anova(fit.1,fit.2)

  if(rocpr==1){
    if(ds){
      print( visreg::visreg(fit.2,"Norm", by="Out", partial=pp, gg=TRUE)  + ggplot2::ylab(latex2exp::TeX('$y_{ij.}$')) + ggplot2::theme_bw() )  # +ggplot2::ylim(0.5,0.71)
    }else{
      print( visreg::visreg(fit.2,"Norm", by="Out", partial=pp, gg=TRUE) +ggplot2::ylim(0.54,0.69) + ggplot2::ylab(latex2exp::TeX('$y_{ij.}$')) + ggplot2::theme_bw() )

    }

  }else{
    print( visreg::visreg(fit.2,"Norm", by="Out", partial=pp, gg=TRUE) +ggplot2::ylim(0.05,0.20)  + ggplot2::ylab(latex2exp::TeX('$y_{ij.}$')) + ggplot2::theme_bw() )
  }

  out <- list()
  out$fit1 <- fit.1
  out$fit2 <- fit.2
  out$aov <- aov_obj
  return(out)

}


IsNormalizingBetter <- function(){

  e <- new.env()
  # ROC values
  data("perf_vals_roc_all", envir=e)
  perfs_all <- e$perf_vals_roc_all
  data("not_normed_perfs_roc_14", envir=e)
  perfs_no <- e$not_normed_perfs_roc_14[, -1]


  pvalues <- matrix(0, nrow=ncol(perfs_no), ncol=4)
  column_names <- c("Mean_SD", "Median_IQR", "Median_MAD", "Min_Max")
  conf_int <- matrix(0, nrow=(ncol(perfs_no)*4), ncol=2)
  mean_vals <- matrix(0, nrow=(ncol(perfs_no)*4), ncol=2)
  normed_better <- matrix(0, nrow=ncol(perfs_no), ncol=4)
  colnames(normed_better) <- column_names

  colnames(pvalues) <- column_names
  methods <-  c()
  for(i in 1:dim(perfs_no)[2]){
    pos <- regexpr("_", colnames(perfs_no)[i]) -1
    method <- substring(colnames(perfs_no)[i], 1, pos)
    methods <- c(methods, method)

    for(j in 1:4){
      k <- (i-1)*4 + j
      test_out <- stats::t.test(perfs_no[ ,i],perfs_all[ ,k] )
      pvalues[i,j] <- test_out$p.value
      conf_int[k,] <- test_out$conf.int[1:2]
      mean_vals[k, ] <- test_out$estimate

      if(test_out$p.value < 0.05){
        normed_better[i, j] <- ifelse((test_out$estimate[2]- test_out$estimate[1]>0), 1, -1)
      }


    }
  }
  methods[which(methods=="FAST")] <- "FAST_ABOD"
  rownames(pvalues) <- methods
  rownames(normed_better) <- methods

  out <- list()
  out$pvalues <- pvalues
  out$confints <- conf_int
  out$mean_vals <- mean_vals
  out$methods <- methods
  out$normed_better <- normed_better

  return(out)
}
