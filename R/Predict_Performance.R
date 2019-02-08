PredictPerformance <- function(ftrs, models){
  ftrs <- as.data.frame(ftrs)
  if( dim(ftrs)[2] > length(models$cols) ){
     x <- ftrs[,models$cols]
  }
  preds <- list()
  pred$cof <- predict(models$cof,newdata=as.matrix(x))
  pred$fabod <- predict(models$fabod,newdata=x)
  pred$inflo <- predict(models$inflo,newdata=x)
  pred$kdeos <- predict(models$kdeos,newdata=x)
  pred$knn <- predict(models$knn,newdata=x)
  pred$knnw <- predict(models$knnw,newdata=x)
  pred$ldf <- predict(models$ldf,newdata=x)
  pred$ldof <- predict(models$ldof,newdata=x)
  pred$lof <- predict(models$lof,newdata=x)
  pred$loof <- predict(models$loof,newdata=x)
  pred$odin <- predict(models$odin,newdata=x)
  pred$simlof <- predict(models$simlof,newdata=x)

  return(pred)
}

#' Train models to predict outlier methods from meta-features of datasets.
#' 
#' @param d If \code{d=1} then we take the Min_Max performance values, if \code{d=2} then performance values from all normalization methods are considered. Input values for  \code{d} are only \code{1,2}.
#' @param p If \code{p=1} then we take binary values based on absolute performance, i.e. if performance \code{ > 0.8}, if \code{p=2} the relative binary performance values are used. Input values for  \code{p} are only \code{1,2}.
#' @param s If \code{s=1} then we train the models on a preferred subset of features. If \code{s=2} the models are trained on all features, which takes considerably longer.  Default value \code{s=1}. Input values for  \code{s} are only \code{1,2}.
#' 
#' @return  The trained randomforest models. 
#'  
#' @examples
#' fit <- TrainModels(1,1,1)
#' 
#' 
TrainModels <- function(d, p, s=1){
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

  if(d==1){
    # ---- ONLY MIN_MAX NORMALIZATION METHOD
    data(features_mm)
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
      data(abs_perfs_mm)
      perfs <- abs_perfs_mm

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_mm)
      perfs <- rel_perfs_0.05_mm

    }

  }else if(d==2){
    # ---- ALL NORMALIZATION METHODS - PERFORMANCE AND FEATURES
    data(features_all)

    if(s==1){
      col_list <- c('SNR', 'OPO_Res_KNOut_95P_1', 'OPO_Out_DenOut_1_3', 'OPO_Den_Out_SD_3', 'OPO_Res_Out_95P_3', 'OPO_LocDenOut_Out_95P_1', 'OPO_GDeg_Out_Mean_1', 'OPO_GComp_PO_Q95_3')
      col_nums <- which(colnames(features_all) %in% col_list )
      ftr_subset <- features_all[ ,col_nums]
    }else if (s==2){
      col_nums <- 1:dim(features_all)[2]
      ftr_subset <- features_all
    }

    if(p==1){
      # absolute performance  1 or 0
      data(abs_perfs_all)
      perfs <- abs_perfs_all

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_all)
      perfs <- rel_perfs_0.05_all
    }
  }
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
  models$cols <- col_nums
  models$d <- d
  models$p <- p
  models$s <- s
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




CrossValidateModels <- function(d, p, s=1, n=5){
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

  if(n > 10){
    stop("Consider n less than or equal to 10.")
  }

  if(d==1){
    # ---- ONLY MIN_MAX NORMALIZATION METHOD
    data(features_mm)
    filenames <- features_mm$filename

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
      data(abs_perfs_mm)
      perfs <- abs_perfs_mm

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_mm)
      perfs <- rel_perfs_0.05_mm

    }

  }else if(d==2){
    # ---- ALL NORMALIZATION METHODS - PERFORMANCE AND FEATURES
    data(features_all)
    filenames <- features_all$filename

    if(s==1){
      col_list <- c('SNR', 'OPO_Res_KNOut_95P_1', 'OPO_Out_DenOut_1_3', 'OPO_Den_Out_SD_3', 'OPO_Res_Out_95P_3', 'OPO_LocDenOut_Out_95P_1', 'OPO_GDeg_Out_Mean_1', 'OPO_GComp_PO_Q95_3')
      col_nums <- which(colnames(features_all) %in% col_list )
      ftr_subset <- features_all[ ,col_nums]

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
      data(abs_perfs_all)
      perfs <- abs_perfs_all

    }else if(p==2){
      # relative performance  1 or 0
      data(rel_perfs_0.05_all)
      perfs <- rel_perfs_0.05_all
    }
  }

  # features are in ftr_subset
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
    for(j in 1:dim(perfs)[2]){
      cat("Fold ", i, " Method " , j, "... \n")
      model <- randomForest::randomForest(trainData, as.factor(trainLabels[ ,j]))
      preds <- predict(model, testData, type="class")
      result_table[i,j] <- sum(preds==testLabels[ ,j])/length(testLabels[ ,j])
      print(paste("Accuracy",  result_table[i,j]) )
    }
  }
  default_accuracy <- apply(perfs, 2, table)*100/dim(perfs)[1]
  default_accuracy <- default_accuracy[1, ]
  out <- list()
  out$def_acc <- default_accuracy
  out$results <-result_table
  out$mean_acc <- apply(result_table, 2, mean)
  return(out)
}
