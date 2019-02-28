#' Cross validates models to predict sensitivity to normalization.
#'
#' @inheritParams SensitivityNormDatasets
#' @param n The number of folds in cross validation.
#'
#' @return A list containing the following results from the cross validated models:
#' \describe{
#'  \item{def_acc}{The default accuracy we get if we predict the method is not good for all instances. This is the percentage of the majority class.}
#'  \item{results}{The \code{n}-fold cross valdation results. }
#'  \item{mean_acc}{The mean \code{n}-fold cross valdation results.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- CrossValidateSensitivityToNorm(1,0.05,10)
#' out$mean_acc
#' }
#'
#' @export

CrossValidateSensitivityToNorm <- function(rocpr=1, xi=0.05, n=10){

  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }

  if(n > 10){
    stop("Consider n less than or equal to 10.")
  }

  e <- new.env()
  if(rocpr==1){
    # ROC values are used
    data(features_all, envir=e)
    filenames <- features_all$filename
    feat <- features_all
    data("perf_vals_roc_all", envir=e)
    perfs <- perf_vals_roc_all

  }else{
    # PR values are used
    data(features_all_pr, envir=e)
    filenames <- features_all_pr$filename
    feat <- features_all_pr
    data("perf_vals_pr_all", envir=e)
    perfs <- perf_vals_pr_all
  }

  col_nums <- 1:dim(feat)[2]
  sds <- apply(feat, 2, sd)
  rm_cols <- c(1, which(sds==0)) # to remove the filename
  col_nums <- setdiff(col_nums,rm_cols)
  ftrs <- feat[ ,-rm_cols]
  cat("Using all", length(col_nums), "features. This will take REALLY LONG!  \n")


  # Compute rangediff
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
  rdiff <-  matrix(0, nrow=dim(perfs)[1], ncol=num_methods)
  rdiff[rangediff > xi] <- 1

  # features are in ftrs
  # performance values in perfs
  result_table <- matrix(0, nrow=n, ncol=dim(perfs)[2])

  # Cross validation on file source as many variants of the same file exist
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

  # Create n equally size folds
  set.seed(1234)
  new_order <- sample(uniq_f_s,length(uniq_f_s))
  folds <- cut(seq(1,length(uniq_f_s)),breaks=n,labels=FALSE)

  cat("Starting", n, "fold cross validation. \n")

  ftr_subset <- apply(ftrs,2,unitize_1)
  # Perform n fold cross validation
  for(i in 1:n){
    testSources <- new_order[which(folds==i,arr.ind=TRUE)]
    testIndices <- which(file_source %in% testSources)

    # Segement your data by fold
    testData <- ftr_subset[testIndices, ]
    trainData <- ftr_subset[-testIndices, ]
    testLabels <- rdiff[testIndices, ]
    trainLabels <- rdiff[-testIndices, ]
    for(j in 1:dim(rdiff)[2]){
      cat("Fold ", i, " Method " , j, "... \n")
      model <- randomForest::randomForest(trainData, as.factor(trainLabels[ ,j]))
      preds <- predict(model, testData, type="class")
      result_table[i,j] <- sum(preds==testLabels[ ,j])/length(testLabels[ ,j])
      print(paste("Accuracy",  result_table[i,j]) )
    }
  }
  default_accuracy <- apply(rdiff, 2, table)*100/dim(rdiff)[1]
  default_accuracy <- default_accuracy[1, ]
  out <- list()
  out$def_acc <- default_accuracy
  out$results <-result_table
  out$mean_acc <- apply(result_table, 2, mean)
  return(out)
}
