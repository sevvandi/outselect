PredictPerformance1 <- function(ftrs, models){
  ftrs <- as.data.frame(ftrs)
  if( dim(ftrs)[2] > length(models$cols) ){
     x <- ftrs[,models$cols]
  }
  predict(models$cof,x)
}

TrainModels <- function(n){
  # n = 1 for absolute performance, i.e. if perf > 0.8 then 1 else 0
  # n = 2 for relative performance with 0.05, i.e. 1 for |perf - max_perf | < 0.05, 0 otherwise
  if((n!=1)|(n!=2)){
    return("Invalid n. n should equal 1 or 2.")
  }

  data(features)
  data(abs_perfs)
  data(rel_perfs_0.05)
  col_list <- c('SNR', 'OPO_Res_KNOut_95P_1', 'OPO_Out_DenOut_1_3', 'OPO_Den_Out_SD_3', 'OPO_Res_Out_95P_3', 'OPO_LocDenOut_Out_95P_1', 'OPO_GDeg_Out_Mean_1', 'OPO_GComp_PO_Q95_3')
  col_nums <- which(colnames(features) %in% col_list )
  ftr_subset <- features[,col_nums]
  if(n==1){
    # Train models
    mod_cof <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,1]) )
    mod_fabod <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,2]) )
    mod_inflo <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,3]) )
    mod_kdeos <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,4]) )
    mod_knn <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,5]) )
    mod_knnw <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,6]) )
    mod_ldf <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,7]) )
    mod_ldof <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,8]) )
    mod_lof <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,9]) )
    mod_loop <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,10]) )
    mod_odin <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,11]) )
    mod_simlof <- randomForest::randomForest(ftr_subset, as.factor(abs_perfs[ ,12]) )
  }else if(n==2){
    # Train models
    mod_cof <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,1]) )
    mod_fabod <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,2]) )
    mod_inflo <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,3]) )
    mod_kdeos <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,4]) )
    mod_knn <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,5]) )
    mod_knnw <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,6]) )
    mod_ldf <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,7]) )
    mod_ldof <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,8]) )
    mod_lof <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,9]) )
    mod_loop <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,10]) )
    mod_odin <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,11]) )
    mod_simlof <- randomForest::randomForest(ftr_subset, as.factor(rel_perfs_0.05[ ,12]) )
  }


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

  return(models)

}

