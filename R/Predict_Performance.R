#' Predicts suitable outlier methods from meta-features of the dataset using trained random forest models.
#' @param ftrs Meta-features of the dataset.
#' @param models The trained random forest models. These are obtained from the function \code{TrainModels}.
#'
#' @return The prediction probabilities for each of the outlier method. For example, a probability of 0.78 for the first outlier method means the probability that the first method is good for this dataset is 0.78.
#'
#' @examples
#' \dontrun{
#' data(Arrhythmia_withoutdupl_05_v05)
#' dat <- Arrhythmia_withoutdupl_05_v05
#' feat <- ComputeMetaFeaturesMM(dat)
#' fit <- TrainModels(1,1,1)
#' out <- PredictPerformance(feat, fit)
#' }

#'
#' @export


PredictPerformance <- function(ftrs, models){
  ftrs <- as.data.frame(ftrs)
  nn <- dim(ftrs)[1]

  d <- models$d
  p <- models$p
  s <- models$s
  col_names <- models$col_names

  x <- ftrs[,col_names]
  if(nn==1){
    # only one record
    if(d==1){
      preds <- matrix(0,nrow=nn, ncol=12)
      colnames(preds) <- c("COF", "FAST_ABOD", "INFLO", "KDEOS", "KNN", "KNNW", "LDF", "LDOF", "LOF", "LOOP", "ODIN", "SIMLOF")

      preds[,1] <- predict(models$cof,newdata=x, type="prob")[2]
      preds[,2]<- predict(models$fabod,newdata=x, type="prob")[2]
      preds[,3] <- predict(models$inflo,newdata=x, type="prob")[2]
      preds[,4] <- predict(models$kdeos,newdata=x, type="prob")[2]
      preds[,5] <- predict(models$knn,newdata=x, type="prob")[2]
      preds[,6] <- predict(models$knnw,newdata=x, type="prob")[2]
      preds[,7] <- predict(models$ldf,newdata=x, type="prob")[2]
      preds[,8] <- predict(models$ldof,newdata=x, type="prob")[2]
      preds[,9]<- predict(models$lof,newdata=x, type="prob")[2]
      preds[,10] <- predict(models$loop,newdata=x, type="prob")[2]
      preds[,11] <- predict(models$odin,newdata=x, type="prob")[2]
      preds[,12] <- predict(models$simlof,newdata=x, type="prob")[2]
    }else{
      # d=2
      preds <- matrix(0,nrow=nn, ncol=8)
      colnames(preds) <- c("Ensemble_Median_IQR", "LOF_Min_Max", "KNN_Median_IQR", "FAST_ABOD_Min_Max", "iForest_Median_IQR",  "KDEOS_Median_IQR", "KDEOS_Min_Max", "LDF_Min_Max")
      preds[,1] <- predict(models$mod_ens_med_iqr,newdata=x, type="prob")[2]
      preds[,2]<- predict(models$mod_lof_min_max,newdata=x, type="prob")[2]
      preds[,3] <- predict(models$mod_knn_med_iqr,newdata=x, type="prob")[2]
      preds[,4] <- predict(models$mod_fabod_min_max,newdata=x, type="prob")[2]
      preds[,5] <- predict(models$mod_iforest_med_iqr,newdata=x, type="prob")[2]
      preds[,6] <- predict(models$mod_kdeos_med_iqr,newdata=x, type="prob")[2]
      preds[,7] <- predict(models$mod_kdeos_min_max,newdata=x, type="prob")[2]
      preds[,8] <- predict(models$mod_ldf_min_max,newdata=x, type="prob")[2]

    }
  }else{
    if(d==1){
      preds <- matrix(0,nrow=nn, ncol=12)
      colnames(preds) <- c("COF", "FAST_ABOD", "INFLO", "KDEOS", "KNN", "KNNW", "LDF", "LDOF", "LOF", "LOOP", "ODIN", "SIMLOF")

      preds[,1] <- predict(models$cof,newdata=x, type="prob")[ ,2]
      preds[,2]<- predict(models$fabod,newdata=x, type="prob")[ ,2]
      preds[,3] <- predict(models$inflo,newdata=x, type="prob")[ ,2]
      preds[,4] <- predict(models$kdeos,newdata=x, type="prob")[ ,2]
      preds[,5] <- predict(models$knn,newdata=x, type="prob")[ ,2]
      preds[,6] <- predict(models$knnw,newdata=x, type="prob")[ ,2]
      preds[,7] <- predict(models$ldf,newdata=x, type="prob")[ ,2]
      preds[,8] <- predict(models$ldof,newdata=x, type="prob")[ ,2]
      preds[,9]<- predict(models$lof,newdata=x, type="prob")[ ,2]
      preds[,10] <- predict(models$loop,newdata=x, type="prob")[ ,2]
      preds[,11] <- predict(models$odin,newdata=x, type="prob")[ ,2]
      preds[,12] <- predict(models$simlof,newdata=x, type="prob")[ ,2]
    }else{
      # d=2
      preds <- matrix(0,nrow=nn, ncol=8)
      colnames(preds) <- c("Ensemble_Median_IQR", "LOF_Min_Max", "KNN_Median_IQR", "FAST_ABOD_Min_Max", "iForest_Median_IQR",  "KDEOS_Median_IQR", "KDEOS_Min_Max", "LDF_Min_Max")
      preds[,1] <- predict(models$mod_ens_med_iqr,newdata=x, type="prob")[ ,2]
      preds[,2]<- predict(models$mod_lof_min_max,newdata=x, type="prob")[ ,2]
      preds[,3] <- predict(models$mod_knn_med_iqr,newdata=x, type="prob")[ ,2]
      preds[,4] <- predict(models$mod_fabod_min_max,newdata=x, type="prob")[ ,2]
      preds[,5] <- predict(models$mod_iforest_med_iqr,newdata=x, type="prob")[ ,2]
      preds[,6] <- predict(models$mod_kdeos_med_iqr,newdata=x, type="prob")[ ,2]
      preds[,7] <- predict(models$mod_kdeos_min_max,newdata=x, type="prob")[ ,2]
      preds[,8] <- predict(models$mod_ldf_min_max,newdata=x, type="prob")[ ,2]

    }
  }



  return(preds)
}

#' Train models to predict outlier methods from meta-features of datasets.
#'
#' @param d If \code{d=1} then we take the Min_Max performance values, if \code{d=2} then performance values from all normalization methods are considered. Input values for  \code{d} are only \code{1,2}.
#' @param p If \code{p=1} then we take binary values based on absolute performance, i.e. if performance \code{ > 0.8}, if \code{p=2} the relative binary performance values are used. Input values for  \code{p} are only \code{1,2}.
#' @param rocpr For all normalization methods, if \code{rocpr=1} then area under \code{ROC} curve is used, if \code{rocpr=2} then area under \code{PR} curve is used. For Min_Max method, i.e. if \code{d=1}, then only \code{ROC} values are used.
#' @param s If \code{s=1} then we train the models on a preferred subset of features. If \code{s=2} the models are trained on all features, which takes considerably longer.  Default value \code{s=1}. Input values for  \code{s} are only \code{1,2}.
#'
#' @return  The trained randomforest models.
#'
#' @examples
#' \dontrun{
#' fit <- TrainModels(1,1,1,1)
#' }
#'
#' @export
#'
TrainModels <- function(d, p, rocpr=1, s=1){
  # d is for meta-data set
  # d = 1 is the min_max performance data set
  # d = 2 is for all normalization methods performance dataset

  # p = 1 for absolute performance, i.e. if perf > 0.8 then 1 else 0
  # p = 2 for relative performance with 0.05, i.e. 1 for |perf - max_perf | < 0.05, 0 otherwise

  # s is for subset of features or all features
  # s = 1 is for a subset of features
  # s = 2 is for all features

  # rocpr is for roc or pr. This is only for all norm methods
  # rocpr = 1 for roc values
  # rocpr = 2 for pr values

  print("Training models on meta-features. This will take some time.")

  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }

  if((p!=1)&(p!=2)){
    stop("Invalid p. p should equal 1 or 2.")
  }

  if((s!=1)&(s!=2)){
    stop("Invalid s. s should equal 1 or 2.")
  }

  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }

  e <- new.env()
  if(d==1){
    # ---- ONLY MIN_MAX NORMALIZATION METHOD
    data(features_mm, envir = e)
    if(s==1){
      col_list <- c('OPO_Res_ResOut_Median_1', 'OPO_Den_Out_95P_1', 'Mean_Entropy_Attr', 'OPO_Res_Out_Mean_1', 'OPO_GDeg_PO_Mean_1',    'IQR_TO_SD_95', 'OPO_GDeg_Out_Mean_1')
      col_nums <- which(colnames(features_mm) %in% col_list )
      ftr_subset <- features_mm[ ,col_nums]
    }else if (s==2){
      # all features
      col_nums <- 1:dim(features_mm)[2]
      ftr_subset <- features_mm
    }

    if(p==1){
      # absolute performance  1 or 0
      data(abs_perfs_mm, envir = e)
      perfs <- abs_perfs_mm

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_mm, envir = e)
      perfs <- rel_perfs_0.05_mm

    }

  }else if(d==2){
    # ---- ALL NORMALIZATION METHODS - PERFORMANCE AND FEATURES
    if(rocpr==1){
      # ROC values
      data(features_all, envir = e)
      data(perf_vals_roc_subset, envir = e)
      perf_vals <- perf_vals_roc_subset
      feat <- features_all
      col_list <- c('OPO_Res_ResOut_Mean_3', 'OPO_GDeg_Out_Mean_3', 'OPO_GComp_PO_Mean_1', 'MEAN_PROD_IQR', 'OPO_Res_Out_Mean_3', 'OPO_Res_ResOut_Max_3', 'OPO_Out_LocDenOut_1_3')
    }else{
      # PR values
      stop("This functionality will be added in the near future.")
      data(features_all_pr, envir = e)
      data(perf_vals_pr_subset, envir = e)
      perf_vals <- perf_vals_pr_subset
      feat <- features_all_pr
      # UPDATE COL_LIST
      col_list <- c(0)
    }


    if(s==1){
      col_nums <- which(colnames(feat) %in% col_list )
      ftr_subset <- feat[ ,col_nums]
    }else if (s==2){
      col_nums <- 1:dim(feat)[2]
      ftr_subset <- feat
    }

    if(p==1){
      # absolute performance  1 or 0
      abs_perf_1_0 <- apply(perf_vals, 2, function(x) ifelse(x>=0.8,1,0))
      perfs <- abs_perf_1_0

    }else if(p==2){
      # relative performance  1 or 0
      relative_perf_eps_0.05 <- EpsilonGood(perf_vals, 0.05)
      perfs <- relative_perf_eps_0.05
    }
  }

  if(d==1){
    # Train models
    mod_cof <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,1]) )
    mod_fabod <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,2]) )
    mod_inflo <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,3]) )
    mod_kdeos <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,4]) )
    mod_knn <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,5]) )
    mod_knnw <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,6]) )
    mod_ldf <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,7]) )
    mod_ldof <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,8]) )
    mod_lof <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,9]) )
    mod_loop <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,10]) )
    mod_odin <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,11]) )
    mod_simlof <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,12]) )

    models <- list()
    models$cof <- mod_cof
    models$fabod <- mod_fabod
    models$inflo <- mod_inflo
    models$kdeos <- mod_kdeos
    models$knn <- mod_knn
    models$knnw <- mod_knnw
    models$ldf <- mod_ldf
    models$ldof <- mod_ldof
    models$lof <- mod_lof
    models$loop <- mod_loop
    models$odin <- mod_odin
    models$simlof <- mod_simlof
  }else{
    # d =2
    mod_ens_med_iqr <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,1]) )
    mod_lof_min_max <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,2]) )
    mod_knn_med_iqr <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,3]) )
    mod_fabod_min_max <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,4]) )
    mod_iforest_med_iqr <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,5]) )
    mod_kdeos_med_iqr <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,6]) )
    mod_kdeos_min_max <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,7]) )
    mod_ldf_min_max <- randomForest::randomForest(ftr_subset, as.factor(perfs[ ,8]) )

    models <- list()
    models$mod_ens_med_iqr <- mod_ens_med_iqr
    models$mod_lof_min_max <- mod_lof_min_max
    models$mod_knn_med_iqr <- mod_knn_med_iqr
    models$mod_fabod_min_max <- mod_fabod_min_max
    models$mod_iforest_med_iqr <- mod_iforest_med_iqr
    models$mod_kdeos_med_iqr <- mod_kdeos_med_iqr
    models$mod_kdeos_min_max <- mod_kdeos_min_max
    models$mod_ldf_min_max <- mod_ldf_min_max
  }

  models$cols <- col_nums
  models$d <- d
  models$p <- p
  models$s <- s
  models$col_names <- col_list
  models$rocpr <- rocpr
  return(models)
}


#' Cross validates models to reproduce results in  journal and conference papers.
#'
#' @inheritParams TrainModels
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
#' out <- CrossValidateModels(1,1,1,1,5)
#' out$mean_acc
#' }
#'
#' @export

CrossValidateModels <- function(d, p, rocpr=1, s=1, n=5){
  # d is for meta-data set
  # d = 1 is the min_max performance data set
  # d = 2 is for all normalization methods performance dataset

  # p = 1 for absolute performance, i.e. if perf > 0.8 then 1 else 0
  # p = 2 for relative performance with 0.05, i.e. 1 for |perf - max_perf | < 0.05, 0 otherwise

  # s is for subset of features or all features
  # s = 1 is for a subset of features
  # s = 2 is for all features

  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }

  if((p!=1)&(p!=2)){
    stop("Invalid p. p should equal 1 or 2.")
  }

  if((s!=1)&(s!=2)){
    stop("Invalid s. s should equal 1 or 2.")
  }

  if((rocpr!=1)&(rocpr!=2)){
    stop("Invalid rocpr. rocpr should equal 1 or 2.")
  }

  if(n > 10){
    stop("Consider n less than or equal to 10.")
  }

  e <- new.env()
  if(d==1){
    # ---- ONLY MIN_MAX NORMALIZATION METHOD
    data(features_mm, envir=e)
    filenames <- features_mm$filename
    data("perf_vals_mm", envir=e)
    perf_best <- apply(perf_vals_mm, 1, which.max)


    if(s==1){
      col_list <- c('OPO_Res_ResOut_Median_1', 'OPO_Den_Out_95P_1', 'Mean_Entropy_Attr', 'OPO_Res_Out_Mean_1', 'OPO_GDeg_PO_Mean_1',    'IQR_TO_SD_95', 'OPO_GDeg_Out_Mean_1')
      col_nums <- which(colnames(features_mm) %in% col_list )
      ftr_subset <- features_mm[ ,col_nums]
    }else if (s==2){
      # all features
      col_nums <- 1:dim(features_mm)[2]
      ftr_subset1 <- features_mm
      sds <- apply(ftr_subset1, 2, sd)
      rm_cols <- c(1, which(sds==0)) # to remove the filename
      col_nums <- setdiff(col_nums,rm_cols)
      ftr_subset <- features_mm[ ,col_nums]
      cat("Using all", length(col_nums), "features. This will take REALLY LONG!  \n")
    }

    if(p==1){
      # absolute performance  1 or 0
      data(abs_perfs_mm, envir=e)
      perfs <- abs_perfs_mm

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_mm, envir=e)
      perfs <- rel_perfs_0.05_mm

    }

  }else if(d==2){
    # ---- ALL NORMALIZATION METHODS - PERFORMANCE AND FEATURES

    if(rocpr==1){
      # ROC values are used
      data(features_all, envir=e)
      filenames <- features_all$filename
      feat <- features_all
      col_list <- c('OPO_Res_ResOut_Mean_3', 'OPO_GDeg_Out_Mean_3', 'OPO_GComp_PO_Mean_1', 'MEAN_PROD_IQR', 'OPO_Res_Out_Mean_3', 'OPO_Res_ResOut_Max_3', 'OPO_Out_LocDenOut_1_3')
      data("perf_vals_roc_subset", envir=e)
      perf_vals <- perf_vals_roc_subset



    }else{
      # PR values are used
      stop("This functionality will be added in the near future.")
      data(features_all_pr, envir=e)
      filenames <- features_all_pr$filename
      feat <- features_all_pr
      data("perf_vals_pr_subset", envir=e)
      perf_vals <- perf_vals_pr_subset
      # NEED TO ADD COL_LIST
      col_list <- c()
    }

    perf_best <- apply(perf_vals, 1, which.max)

    if(s==1){
      col_nums <- which(colnames(features_all) %in% col_list )
      ftr_subset <- feat[ ,col_nums]

    }else if (s==2){
      col_nums <- 1:dim(features_all)[2]
      ftr_subset1 <- features_all
      sds <- apply(ftr_subset1, 2, sd)
      rm_cols <- c(1, which(sds==0)) # to remove the filename
      col_nums <- setdiff(col_nums,rm_cols)
      ftr_subset <- features_all[ ,-rm_cols]
      cat("Using all", length(col_nums), "features. This will take REALLY LONG!  \n")
    }

    if(p==1){
      # absolute performance  1 or 0
      abs_perf_1_0 <- apply(perf_vals, 2, function(x) ifelse(x>=0.8,1,0))
      perfs <- abs_perf_1_0

    }else if(p==2){
      # relative performance  1 or 0
      relative_perf_eps_0.05 <- EpsilonGood(perf_vals, 0.05)
      perfs <- relative_perf_eps_0.05
    }

  }

  # features are in ftr_subset
  # performance values in perfs
  result_table <- matrix(0, nrow=n, ncol=dim(perfs)[2])
  pred_algs <- matrix(0, nrow=length(filenames), ncol=dim(perfs)[2])
  colnames(pred_algs) <- colnames(perfs)
  pred_best <- matrix(0, nrow=length(filenames), ncol=dim(perfs)[2])
  colnames(pred_best) <- colnames(perfs)

  # Cross validation on file source as many variants of the same file exist
  file_source <-GetFileSources(filenames)

  uniq_f_s <- unique(file_source)

  # Create n equally size folds
  set.seed(1234)
  new_order <- sample(uniq_f_s,length(uniq_f_s))
  folds <- cut(seq(1,length(uniq_f_s)),breaks=n,labels=FALSE)

  cat("Starting", n, "fold cross validation. \n")

  ftr_subset <- apply(ftr_subset,2,unitize_1)
  # Perform n fold cross validation
  for(i in 1:n){
    testSources <- new_order[which(folds==i,arr.ind=TRUE)]
    testIndices <- which(file_source %in% testSources)

    # Segement your data by fold
    testData <- ftr_subset[testIndices, ]
    trainData <- ftr_subset[-testIndices, ]
    testLabels <- perfs[testIndices, ]
    trainLabels <- perfs[-testIndices, ]

    # Train model for best method
    testLabBest <- perf_best[testIndices]
    trainLabBest <- perf_best[-testIndices]
    modelBest <-randomForest::randomForest(trainData, as.factor(paste(trainLabBest)))
    pred_best[testIndices,] <- predict(modelBest, testData, type="prob" )

    for(j in 1:dim(perfs)[2]){
      cat("Fold ", i, " Method " , j, "... \n")
      model <- randomForest::randomForest(trainData, as.factor(trainLabels[ ,j]))
      preds <- predict(model, testData, type="class")
      pred_algs[testIndices,j] <- predict(model, testData, type="prob")[ ,2]
      result_table[i,j] <- sum(preds==testLabels[ ,j])/length(testLabels[ ,j])
      print(paste("Accuracy",  result_table[i,j]) )
    }
  }
  default_accuracy <- apply(perfs, 2, table)*100/dim(perfs)[1]
  default_accuracy <- apply(default_accuracy,2, max)
  out <- list()
  out$def_acc <- default_accuracy
  out$results <-result_table
  out$mean_acc <- apply(result_table, 2, mean)
  out$algo_preds <- pred_algs
  out$pred_best <- pred_best
  return(out)
}



#' Cross validates SVM with instance space coordinates.
#'
#'@inheritParams CrossValidateModels
#'
#'@return A list with the following components:
#'\describe{
#'  \item{def_acc}{The default accuracy we get if we predict the method is not good for all instances. This is the percentage of the majority class.}
#'  \item{results}{The \code{n}-fold cross valdation results. }
#'  \item{mean_acc}{The mean \code{n}-fold cross valdation results.}
#'  \item{d}{Which instance space is cross validated. If \code{d=1}, then it is MIN_MAX normalization instance space, if \code{d=2} then all normalization methods are used. }
#'  \item{coordinates}{The instance space coordinates.}
#' }
#'
#'@examples
#'\dontrun{ out <- CrossValidateSVM(1,5) }
#'
#'@export
#'
CrossValidateSVM <- function(d=1,n=5){
  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }

  if(n > 10){
    stop("Consider n less than or equal to 10.")
  }
  e <- new.env()

  if(d==1){
    # coordinates for min-max normalization
    data("dat_4_svm_mm", envir=e)
    coordinates <- dat_4_svm_mm[, c(2,3)]
    perfs <- dat_4_svm_mm[ ,4:15]
    filenames <- dat_4_svm_mm[, 1]
    cst <- 150
    gmv <- 0.75
  }else{
    data("data_4_svm_mix_5K", envir=e)
    coordinates <- data_4_svm_mix_5K[, c(2,3)]
    perfs <- data_4_svm_mix_5K[ ,4:11]
    filenames <- data_4_svm_mix_5K[, 1]
    cst <- 50
    gmv <- 0.25

  }

  # Cross validation on file source as many variants of the same file exist
  file_source <-GetFileSources(filenames)
  uniq_f_s <- unique(file_source)

  # Create n equally size folds
  set.seed(1234)
  new_order <- sample(uniq_f_s,length(uniq_f_s))
  folds <- cut(seq(1,length(uniq_f_s)),breaks=n,labels=FALSE)

  cat("Starting", n, "fold cross validation. \n")

  result_table <- matrix(0, nrow=n, ncol=dim(perfs)[2])

  # Perform n fold cross validation
  for(i in 1:n){
    testSources <- new_order[which(folds==i,arr.ind=TRUE)]
    testIndices <- which(file_source %in% testSources)

    # Segement your data by fold
    testData <- coordinates[testIndices, ]
    trainData <- coordinates[-testIndices, ]
    testLabels <- perfs[testIndices, ]
    trainLabels <- perfs[-testIndices, ]
    for(j in 1:dim(perfs)[2]){
      cat("Fold ", i, " Method " , j, "... \n")
      model <- e1071::svm(trainData, as.factor(trainLabels[ ,j]), cost = cst, gamma=gmv)
      preds <- predict(model, testData, type="class")
      result_table[i,j] <- sum(preds==testLabels[ ,j])/length(testLabels[ ,j])
      print(paste("Accuracy",  result_table[i,j]) )
    }
  }
  default_accuracy <- apply(perfs, 2, table)*100/dim(perfs)[1]
  default_accuracy <- apply(default_accuracy,2, max)
  out <- list()
  out$def_acc <- default_accuracy
  out$results <-result_table
  out$mean_acc <- apply(result_table, 2, mean)*100
  out$d <- d
  out$coordinates <- coordinates
  return(out)
}

#' Trains an SVM on instance space coordinates and draws the instance space.
#'
#' @inheritParams CrossValidateSVM
#' @param vis If \code{TRUE} then the instance space is plotted.
#'
#' @return A list with the following components:
#'\describe{
#'  \item{preds10}{The predictions as \code{1} or \code{0} for each instance and each algorithm.}
#'  \item{predprobs}{The prediction probabilities for each each instance and each algorithm.}
#'  \item{algorithms}{The algorithm best suited for each instance. For some instances, no algorithms work well.}
#'  \item{d}{Which instance space is cross validated. If \code{d=1}, then it is MIN_MAX normalization instance space, if \code{d=2} then all normalization methods are used. }
#'  \item{coordinates}{The instance space coordinates.}
#' }
#'
#'@examples
#'\dontrun{ svmout <- InstSpace(d=1)}
#'
#'@export
#'

InstSpace <- function(d=1, vis=FALSE){
  if((d!=1)&(d!=2)){
    stop("Invalid d. d should equal 1 or 2.")
  }
  e <- new.env()
  if(d==1){
    # Instance space using MIN_MAX
    data("dat_4_svm_mm", envir=e)
    filenames <- dat_4_svm_mm[ ,1]
    xx <- dat_4_svm_mm[, 2:3]
    perfs <- dat_4_svm_mm[ , 4:15]
    cst <- 150
    gmv <- 0.75
    num_methods <- dim(perfs)[2]
    print("Training 12 SVMS for outlier detection methods. This will take some time.")



  }else{
    # Instance space using all normalization methods
    data("data_4_svm_mix_5K", envir=e)
    filenames <- data_4_svm_mix_5K[ ,1]
    xx <- data_4_svm_mix_5K[, 2:3]
    perfs <- data_4_svm_mix_5K[ , 4:11]
    cst <- 50
    num_methods <- dim(perfs)[2]
    gmv <-  0.25 # 1/num_methods #

    print("Training 8 SVMS for outlier detection methods. This will take some time.")  }

  # Train an SVM
  preds.all <- matrix(0,nrow=nrow(xx), ncol=num_methods)
  preds.all.1.0 <- matrix(0,nrow=nrow(xx), ncol=num_methods)
  for(kk in 1:num_methods){
    dd2 <- cbind.data.frame(xx, as.factor(perfs[,kk]))
    colnames(dd2) <- c("x", "y","Label" )

    # Class weights - BEGIN
    costs <- table(dd2$Label)
    costs[1] <- 1
    costs[2] <- c(1,round(table(dd2$Label)[1]/table(dd2$Label)[2]))[1]
    # Class weights - END

    mod.big <- e1071::svm(Label~.,  data=dd2, probability=TRUE, class.weights=costs, scale=FALSE, cost=cst, gamma=gmv )
    pf.all <- predict(mod.big, dd2, probability=TRUE)
    col <- which(colnames(attr(pf.all, "probabilities"))=="1")
    preds.all[,kk] <- attr(pf.all, "probabilities")[,col]
    preds.all.1.0[,kk] <- ifelse(preds.all[,kk] >=0.5,1,0)
  }
  colnames(preds.all.1.0)<-  colnames(perfs)


  tt3 <- apply(preds.all.1.0,1,sum)
  multiples <- which(tt3>1)
  pp.multi <- apply(preds.all[multiples,],1,function(x) which.max(x))

  pred.m <- rep(0,dim(preds.all.1.0)[1])
  pred.m[multiples] <-pp.multi
  qq <- apply(preds.all.1.0[-multiples,],1,function(x) ifelse(sum(x)==1,which(x==1), 0 ) )
  pred.m[is.na(pred.m)] <- 0

  pred.m[-multiples] <- qq
  head(pred.m)
  algorithms <- pred.m
  if(d==1){
    algorithms[pred.m==0] <- "None"
    algorithms[pred.m==1] <- "COF"
    algorithms[pred.m==2] <- "FAST ABOD"
    algorithms[pred.m==3] <- "INFLO" # "None" #
    algorithms[pred.m==4] <- "KDEOS"
    algorithms[pred.m==5] <- "KNN"
    algorithms[pred.m==6] <- "KNNW"
    algorithms[pred.m==7] <-  "LDF" # "None" #
    algorithms[pred.m==8] <-  "LDOF" #
    algorithms[pred.m==9] <-  "LOF" #
    algorithms[pred.m==10] <-  "LOOP" #
    algorithms[pred.m==11] <-  "ODIN" #
    algorithms[pred.m==12] <-  "SIMLOF" #
  }else{
    algorithms[pred.m==0] <- "None"
    algorithms[pred.m==1] <- "Ensemble_Median_IQR"
    algorithms[pred.m==2] <- "LOF_Min_Max"
    algorithms[pred.m==3] <- "KNN_Median_IQR"
    algorithms[pred.m==4] <- "FAST_ABOD_Min_Max"
    algorithms[pred.m==5] <- "iForest_Median_IQR"
    algorithms[pred.m==6] <- "KDEOS_Median_IQR"
    algorithms[pred.m==7] <-  "KDEOS_Min_Max"
    algorithms[pred.m==8] <-  "LDF_Min_Max"
  }


  if(vis){
    algorithms <- as.factor(algorithms)
    print( ggplot2::ggplot(data=xx, ggplot2::aes(x,y))+ ggplot2::geom_point(ggplot2::aes(color=algorithms, shape=algorithms)) + ggplot2::scale_shape_manual(values=1:nlevels(algorithms)) +  ggplot2::theme_bw() )
  }

  out <- list()
  out$preds10 <- preds.all.1.0
  out$predprobs <- preds.all
  out$algorithms <- algorithms
  out$filenames <- filenames
  out$d <- d
  out$coordinates <- xx
  return(out)
}

#' Plots a new instance in the instance space.
#'
#' @param svm_out The output of the trained svm using function \code{InstSpace}
#' @param feat The features of the new instance/dataset. This can be computed using \code{ComputeMetaFeaturesAll} or \code{ComputeMetaFeaturesMM}.
#' @param vis If \code{TRUE} then the instance space along with the new instance is plotted.
#'
#' @return new_coords Coodinates of the new instance in the instance space
#'
#' @examples
#' \dontrun{
#' data(Arrhythmia_withoutdupl_05_v05)
#' dat <- Arrhythmia_withoutdupl_05_v05
#' feat <- ComputeMetaFeaturesMM(dat)
#' svmout <- InstSpace(d=1)
#' PlotNewInstance(svmout, feat, vis=FALSE)
#' }
#'
#'@export
PlotNewInstance <- function(svm_out, feat, vis=TRUE){
# ---- NEED TO UPDATE THIS FUNCTION
# ---- FEAT NEED TO BE TRANSFORMED FIRST
  d <- svm_out$d
  coordinates <- svm_out$coordinates
  algorithms <- svm_out$algorithms
  feat <- as.data.frame(feat)

  if(d==1){
    # MIN_MAX NORMALIZATION
    proj_mat <- matrix( c(-0.0460, 0.1202, -0.0862, -0.0938, 0.1854, 0.1737, 0.3543, -0.2847, 0.0378, -0.2078,-0.2025, -0.0822, 0.1845, -0.1325), nrow=2 )
    col_names <- c("OPO_Res_ResOut_Median_1", "OPO_Den_Out_95P_1", "Mean_Entropy_Attr", "OPO_Res_Out_Mean_1", "OPO_GDeg_PO_Mean_1", "IQR_TO_SD_95", "OPO_GDeg_Out_Mean_1")
    x_bef <- feat[ , col_names]
   # xx <- renormalize(x_bef, d)
    xx <- x_bef
  }else{
    # ALL NORMALIZATION METHODS
    col_names <- c("OPO_Res_ResOut_Mean_3","OPO_GDeg_Out_Mean_3", "OPO_GComp_PO_Mean_1", "MEAN_PROD_IQR", "OPO_Res_Out_Mean_3", "OPO_Res_ResOut_Max_3",   "OPO_Out_LocDenOut_1_3" )
    x_bef <- feat[ , col_names]
    xx <- renormalize(x_bef, d)
    proj_mat <- matrix( c(-0.00292627076652680,-0.632192408787632,-0.569663704435377,-0.174829830230183,0.0268201438117375,-0.290501129270193,-0.139810891656474, 0.332268685436933,-0.224517922081636,-0.0629031205759270,0.182180599272572,0.543197310499069,0.427965822151968,0.389217122011117 ), nrow=2, byrow=TRUE )

  }
  new_coords <- t(proj_mat%*%t(xx))
  colnames(new_coords) <- colnames(coordinates)

  if(vis){
    algorithms <- as.factor(algorithms)
    print(ggplot2::ggplot(data=coordinates, ggplot2::aes(x,y))+ ggplot2::geom_point(ggplot2::aes(color=algorithms)) + ggplot2::geom_point(ggplot2::aes(x=new_coords[1],y=new_coords[2]), color="black", shape=17, size=4)  + ggplot2::theme_bw() )
  }

  return(new_coords)
}





