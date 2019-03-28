#' Clusters outlier detection methods according to the correlation matrix for min-max and median-iqr normalization methods, using the R package \code{cluster}.
#'
#' @param kk The number of clusters
#' @param rocpr If \code{rocpr = 1}, then area under ROC cure is used, if \code{rocpr = 2} then area under PR curve is used.
#' @param vis If \code{TRUE} then the resulting graphs are plotted.
#'
#' @return A \code{pam.object} described in the R package \code{cluster}.
#' @seealso  \code{\link[cluster]{pam.object}}.
#'
#' @examples
#' cls <- ClusterMethods(8)
#' which(cls$clustering==2)
#' which(cls$clustering==4)
#' which(cls$clustering==5)
#' which(cls$clustering==6)
#'
#'@importFrom graphics plot
#'@export
ClusterMethods <- function( kk, rocpr=1, vis=FALSE){
  # kk number of clusters

  # rocpr = 1 for roc values, rocpr =2 for pr values

  if(kk <2){
    stop("kk needs to be greater than 1.")
  }
  if((rocpr!=1)&(rocpr!=2)){
    stop("rocpr can only be 1 or 2.")
  }
  e <- new.env()
  if(rocpr==1){
    data(perf_vals_roc_all, envir=e)
    dat <- perf_vals_roc_all[ ,2*1:28]
  }else{
    data(perf_vals_pr_all, envir=e)
    dat <- perf_vals_pr_all[ ,2*1:28]
  }
  cor_obj <- stats::cor(dat)
  dissimilarity <- 1 - cor_obj
  distance <- stats::as.dist(dissimilarity)
  pam_obj <- cluster::pam(distance, k=kk)
  if(vis){
    print(plot(pam(distance, k=kk)))
  }
  return(pam_obj)
}

#' Computes correlation between AUCROC performance values and AUCPR performance values
#'
#' @param m If \code{1}, then \code{pearson}, if \code{2}, then \code{kendall} or if \code{3}, then \code{spearman} correlation is performed. Default value is \code{1}.
#'
#' @return  The correlation values for each performance column.
#'
#' @examples
#' \dontrun{
#' corobj <- CorrPrRoc()
#' }
#'
#' @export
CorrPrRoc <- function(m=1){
  if(m==1){
    meth ="pearson"
  }else if(m==2){
    meth ="kendall"
  }else if(m==3){
    meth ="spearman"
  }else{
    stop("Invalid value for m. m needs to be 1, 2 or 3.")
  }

  e <- new.env()
  data("features_all", envir=e)
  data("features_all_pr",envir=e)
  fnames_roc <- features_all[,1]
  fnames_pr <- features_all_pr[,1]
  roc_list <- which(fnames_roc %in% fnames_pr)
  pr_list <- which(fnames_pr %in% fnames_roc)
  data("perf_vals_pr_all", envir=e)
  data("perf_vals_roc_all", envir=e)
  cor_vals <- rep(0,56)
  for(i in 1:56){
    if(identical(colnames(perf_vals_roc_all)[i], colnames(perf_vals_pr_all)[i] ) ){
      cor_vals[i] <- stats::cor(perf_vals_roc_all[roc_list,i], perf_vals_pr_all[pr_list,i], use="complete.obs", method=meth)
    }
  }
  return(cor_vals)
}

#'@export
DifficultyDiversitySpace <- function(d=2, rocpr=1){
  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }
  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }
  e <- new.env()
  if(d==1){
    # MIN-MAX normalization
    data("perf_vals_mm", envir=e)
    perfs <- perf_vals_mm
    data("filenames_mm", envir=e)
    filenames <- filenames_mm
  }else{
    # All normalization methods
    if(rocpr==1){
      # ROC values
      data("perf_vals_roc_subset", envir=e)
      perfs <- perf_vals_roc_subset
      data("filenames_roc", envir=e)
      filenames <- filenames_roc
    }else{
      # PR values
      data("perf_vals_pr_subset", envir=e)
      perfs <- perf_vals_pr_subset
      data("filenames_pr", envir=e)
      filenames <- filenames_pr
    }
  }

  # difficulty is 1 - average performance
  mean_vals <- 1- apply(perfs, 1, mean)
  sd_vals <- apply(perfs, 1, sd)

  file_source <- GetFileSources(filenames)
  df <- cbind.data.frame(filenames, file_source, mean_vals, sd_vals)
  colnames(df) <- c("filename", "source", "difficulty", "diversity" )
  source_avg <- stats::aggregate(df$difficulty, by=list(df$source), FUN=mean)
  source_std <- stats::aggregate(df$diversity, by=list(df$source), FUN=mean)

  source_df <- cbind.data.frame(source_avg, source_std[ ,2])
  colnames(source_df) <- c("source", "source_difficulty", "source_diversity")

  out <- list()
  out$dfall <- df
  out$dfsrc <- source_df
  return(out)
}


GetSourceRepository <- function(d=2){
  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }
  e <- new.env()
  if(d==1){
    # Min-Max normalization
    data("filenames_mm", envir=e)
    filenames <- e$filenames_mm
  }else{
    # Mix of norm and outlier methods
    data("filenames_roc", envir=e)
    filenames <- e$filenames_roc
  }
  wierd.list <- c()
  file_repo_src <-c()
  # Campos source = 1
  # Goldstein source = 2
  # UCI source = 3
  for(i in 1:length(filenames)){
    fname <- filenames[i]
    regobj1 <- regexpr("_C", fname)  # UCI source = 3
    regobj2 <- regexpr("_withoutdupl", fname)  # Campos source = 1
    regobj3 <- regexpr("unsupervised", fname)  # Goldstein source = 2

    if(regobj2[1]>0){
      src <- 1
    }else if(regobj1[1]>0){
      src <- 3
    }else if(regobj3[1]>0){
      src <- 2
    }else{
      wierd.list <- c(wierd.list, i)
      src <- 0
    }
    file_repo_src <- c(file_repo_src, src)
  }
  file_repo_details <- cbind.data.frame(filenames, file_repo_src)
  return(file_repo_details)
}
