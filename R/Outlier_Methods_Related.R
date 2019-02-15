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
ClusterMethods <- function( kk, rocpr=1, vis=FALSE){
  # kk number of clusters

  # rocpr = 1 for roc values, rocpr =2 for pr values

  if(kk <2){
    stop("kk needs to be greater than 1.")
  }
  if(rocpr==1){
    data(perf_vals_roc_all)
    dat <- perf_vals_roc_all[ ,2*1:28]
  }else{
    stop("For later - for PR values")
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
