
library(MASS)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(caret)
library(nnet)
library(pROC)

library(foreach)
library(doParallel)

source('simulations/00_sim_funs.R')
source('simulations/00_lmm_funs.R')
#---- CLARITY ----
library(Clarity)




t <- 10
d <- 30

n <-  500 # 200
nU <- 2
nC <- 4
nX <- 100 # 200
nbins <- 100
clusters = 7
rho_time = 0.8
sigma_x = 0.0375
sigma_c = 0.12
mean_effect = 0.25

n_cores <- 4
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

#------------------- Mean sensitivity analysis -------------------
means <- seq(from = 0, to = 1, length.out = nbins )


aucs_mean <- foreach(i = 1:t, .combine = rbind) %dopar% {
  library(dplyr)
  results <- tibble::tibble(
    Seed = numeric(),
    Dimension = numeric(),
    Method = character(),
    AUC = numeric(),
    Mean = numeric()
  )
  for (m in means){
  
    set.seed(i)
    out <- sim_data(
      n = n, nU = nU, nC = nC, nX = nX, clusters = clusters, rho_time = rho_time,
      sigma_x = sigma_x, sigma_c = sigma_c, mean_effect = m
    )
    
    X1 <- out$X1  
    X2 <- out$X2
    C <- out$UC[,paste0('C',1:nC)]
    anomalies <- out$anomalies
    
    all_data <- rbind(X1, X2)
    
    ### lmm
    lmm_input <- cbind(all_data, 'Timepoint' = rep(c(1,2), each = n), rbind(C,C), studyid = rep(1:n, times = 2))
    lmm <- metab.mixed.model(as.data.frame(lmm_input), columns_with_metabolite_data = paste0('X', 1:nX), covariables = paste0('C',1:3))
    lmm_p_vals <- lmm$P_mixed
    lmm_p_vals <- p.adjust(lmm_p_vals, method = c('BH'))
    
    ### Clarity 
    scaled <- scale(all_data)
    
    X1_cov <- cov(scaled[1:n,])
    X2_cov <- cov(scaled[(n +1):(2*n),])
    
    X1_cov <- Clarity::c_fixdiagonals(X1_cov)
    X2_cov <- Clarity::c_fixdiagonals(X2_cov)
    X1_scan_cov  <- Clarity::Clarity_Scan(X1_cov, kmax = d,verbose=FALSE)
    ## Rotate the baseline to best-fit the endpoint
    X2_pred_cov <- Clarity::Clarity_Predict(X2_cov, X1_scan_cov)
    
    pers_cov <- Clarity::Clarity_Persistence(X2_pred_cov)
    
    ## Combine LMM and Clarity p-vals in multinomial logistic regression
    pvaldf=data.frame(lmm_p_vals=lmm_p_vals, clarity_p_vals=pers_cov[,10],anomalies=anomalies)
    pvaldf$anomalies <- as.factor(as.character(as.numeric(pvaldf$anomalies)))
    pvaldf$anomalies <- relevel(pvaldf$anomalies, ref = "0")
    sapply(1:100,function(rep){
      set.seed(rep)
      index <- caret::createDataPartition(pvaldf$anomalies, p = .70, list = FALSE)
      train <- pvaldf[index,]
      test <- pvaldf[-index,]
      
      lmm_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals, data = train)
      lmm_multinom_pred <- predict(lmm_multinom_model, newdata = test, "class", type = "prob")
      auc_lmm <- as.numeric(pROC::multiclass.roc(test$anomalies, lmm_multinom_pred)$auc)
      clarity_multinom_model <- nnet::multinom(anomalies ~ clarity_p_vals, data = train)
      clarity_multinom_pred <- predict(clarity_multinom_model, newdata = test, "class", type = "prob")
      auc_clarity <- as.numeric(pROC::multiclass.roc(test$anomalies, clarity_multinom_pred)$auc)
      comb_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals + clarity_p_vals, data = train)
      comb_multinom_pred <- predict(comb_multinom_model, newdata = test, "class", type = "prob")
      auc_comb <- as.numeric(pROC::multiclass.roc(test$anomalies, comb_multinom_pred)$auc)
      return(c(auc_lmm, auc_clarity, auc_comb))
    }) -> aucs_matrix
    auc_lmm <- mean(aucs_matrix[1,])
    auc_clarity <- mean(aucs_matrix[2,])
    auc_comb <- mean(aucs_matrix[3,])
    
    results <- results %>% add_row(
      Seed = t, 
      Dimension = NA, 
      Method = 'LMM',
      AUC = auc_lmm,
      Mean = m
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Clarity',
      AUC = auc_clarity,
      Mean = m
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Combined',
      AUC = auc_comb,
      Mean = m
    )    
    
  }
  return(results)
}
   
write.csv(aucs_mean, 'simulations/results/mean_multiway_aucs.csv')

############# Rho time sensitivity analysis ################
rho_times <- seq(from = 0, to = 1, length.out = nbins)


aucs_rho <- foreach(i = 1:t, .combine = rbind) %dopar% {
  library(dplyr)
  results <- tibble::tibble(
    Seed = numeric(),
    Dimension = numeric(),
    Method = character(),
    AUC = numeric(),
    Rho = numeric(),
  )
  for (r in rho_times){
    set.seed(i)
    out <- sim_data(
      n = n, nU = nU, nC = nC, nX = nX, clusters = clusters, rho_time = r,
      sigma_x = sigma_x, sigma_c = sigma_c, mean_effect = mean_effect
    )
    
    X1 <- out$X1  
    X2 <- out$X2
    C <- out$UC[,paste0('C',1:nC)]
    anomalies <- out$anomalies
    
#    anomalies  <- as.logical(anomalies)
    
    all_data <- rbind(X1, X2)
    
    ### lmm
    lmm_input <- cbind(all_data, 'Timepoint' = rep(c(1,2), each = n), rbind(C,C), studyid = rep(1:n, times = 2))
    lmm <- metab.mixed.model(as.data.frame(lmm_input), columns_with_metabolite_data = paste0('X', 1:nX), covariables = paste0('C',1:3))
    lmm_p_vals <- lmm$P_mixed
    lmm_p_vals <- p.adjust(lmm_p_vals, method = c('BH'))
    
    ### Clarity 
    scaled <- scale(all_data)
    
    X1_cov <- cov(scaled[1:n,])
    X2_cov <- cov(scaled[(n +1):(2*n),])
    
    X1_cov <- Clarity::c_fixdiagonals(X1_cov)
    X2_cov <- Clarity::c_fixdiagonals(X2_cov)
    X1_scan_cov  <- Clarity::Clarity_Scan(X1_cov, kmax = d,verbose=FALSE)
    ## Rotate the baseline to best-fit the endpoint
    X2_pred_cov <- Clarity::Clarity_Predict(X2_cov, X1_scan_cov)
    
    pers_cov <- Clarity::Clarity_Persistence(X2_pred_cov)
  
    ## Combine LMM and Clarity p-vals in multinomial logistic regression
    pvaldf=data.frame(lmm_p_vals=lmm_p_vals, clarity_p_vals=pers_cov[,10],anomalies=anomalies)
    pvaldf$anomalies <- as.factor(as.character(as.numeric(pvaldf$anomalies)))
    pvaldf$anomalies <- relevel(pvaldf$anomalies, ref = "0")
    sapply(1:100,function(rep){
      set.seed(rep)
      index <- caret::createDataPartition(pvaldf$anomalies, p = .70, list = FALSE)
      train <- pvaldf[index,]
      test <- pvaldf[-index,]
      
      lmm_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals, data = train)
      lmm_multinom_pred <- predict(lmm_multinom_model, newdata = test, "class", type = "prob")
      auc_lmm <- as.numeric(pROC::multiclass.roc(test$anomalies, lmm_multinom_pred)$auc)
      clarity_multinom_model <- nnet::multinom(anomalies ~ clarity_p_vals, data = train)
      clarity_multinom_pred <- predict(clarity_multinom_model, newdata = test, "class", type = "prob")
      auc_clarity <- as.numeric(pROC::multiclass.roc(test$anomalies, clarity_multinom_pred)$auc)
      comb_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals + clarity_p_vals, data = train)
      comb_multinom_pred <- predict(comb_multinom_model, newdata = test, "class", type = "prob")
      auc_comb <- as.numeric(pROC::multiclass.roc(test$anomalies, comb_multinom_pred)$auc)
      return(c(auc_lmm, auc_clarity, auc_comb))
    }) -> aucs_matrix
    auc_lmm <- mean(aucs_matrix[1,])
    auc_clarity <- mean(aucs_matrix[2,])
    auc_comb <- mean(aucs_matrix[3,])
    
    results <- results %>% add_row(
      Seed = t, 
      Dimension = NA, 
      Method = 'LMM',
      AUC = auc_lmm,
      Rho = r
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Clarity',
      AUC = auc_clarity,
      Rho = r
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Combined',
      AUC = auc_comb,
      Rho = r
    )    
    
    
  }
  return(results)
}

write.csv(aucs_rho, 'simulations/results/rho_multiway_aucs.csv')



############# Cluster sensitivity analysis ################
clusters <- as.numeric(5:20)


aucs_cluster <- foreach(i = 1:t, .combine = rbind) %dopar% {
  library(dplyr)
  results <- tibble::tibble(
    Seed = numeric(),
    Dimension = numeric(),
    Method = character(),
    AUC = numeric(),
    Cluster = numeric(),
  )
  for (c in clusters){
    set.seed(i)
    out <- sim_data(
      n = n, nU = nU, nC = nC, nX = nX, clusters = c, rho_time = rho_time,
      sigma_x = sigma_x, sigma_c = sigma_c, mean_effect = mean_effect
    )
    
    X1 <- out$X1  
    X2 <- out$X2
    C <- out$UC[,paste0('C',1:nC)]
    anomalies <- out$anomalies
    
    all_data <- rbind(X1, X2)
    
    ### lmm
    lmm_input <- cbind(all_data, 'Timepoint' = rep(c(1,2), each = n), rbind(C,C), studyid = rep(1:n, times = 2))
    lmm <- metab.mixed.model(as.data.frame(lmm_input), columns_with_metabolite_data = paste0('X', 1:nX), covariables = paste0('C',1:3))
    lmm_p_vals <- lmm$P_mixed
    lmm_p_vals <- p.adjust(lmm_p_vals, method = c('BH'))
    
    ### Clarity 
    scaled <- scale(all_data)
    
    X1_cov <- cov(scaled[1:n,])
    X2_cov <- cov(scaled[(n +1):(2*n),])
    
    X1_cov <- Clarity::c_fixdiagonals(X1_cov)
    X2_cov <- Clarity::c_fixdiagonals(X2_cov)
    X1_scan_cov  <- Clarity::Clarity_Scan(X1_cov, kmax = d,verbose=FALSE)
    ## Rotate the baseline to best-fit the endpoint
    X2_pred_cov <- Clarity::Clarity_Predict(X2_cov, X1_scan_cov)
    
    pers_cov <- Clarity::Clarity_Persistence(X2_pred_cov)
    
## Combine LMM and Clarity p-vals in multinomial logistic regression
    pvaldf=data.frame(lmm_p_vals=lmm_p_vals, clarity_p_vals=pers_cov[,10],anomalies=anomalies)
    pvaldf$anomalies <- as.factor(as.character(as.numeric(pvaldf$anomalies)))
    pvaldf$anomalies <- relevel(pvaldf$anomalies, ref = "0")
    sapply(1:100,function(rep){
      set.seed(rep)
      index <- caret::createDataPartition(pvaldf$anomalies, p = .70, list = FALSE)
      train <- pvaldf[index,]
      test <- pvaldf[-index,]
    
      lmm_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals, data = train)
      lmm_multinom_pred <- predict(lmm_multinom_model, newdata = test, "class", type = "prob")
      auc_lmm <- as.numeric(pROC::multiclass.roc(test$anomalies, lmm_multinom_pred)$auc)
      clarity_multinom_model <- nnet::multinom(anomalies ~ clarity_p_vals, data = train)
      clarity_multinom_pred <- predict(clarity_multinom_model, newdata = test, "class", type = "prob")
      auc_clarity <- as.numeric(pROC::multiclass.roc(test$anomalies, clarity_multinom_pred)$auc)
      comb_multinom_model <- nnet::multinom(anomalies ~ lmm_p_vals + clarity_p_vals, data = train)
      comb_multinom_pred <- predict(comb_multinom_model, newdata = test, "class", type = "prob")
      auc_comb <- as.numeric(pROC::multiclass.roc(test$anomalies, comb_multinom_pred)$auc)
      return(c(auc_lmm, auc_clarity, auc_comb))
    }) -> aucs_matrix
    auc_lmm <- mean(aucs_matrix[1,])
    auc_clarity <- mean(aucs_matrix[2,])
    auc_comb <- mean(aucs_matrix[3,])
    
    results <- results %>% add_row(
      Seed = t, 
      Dimension = NA, 
      Method = 'LMM',
      AUC = auc_lmm,
      Cluster = c
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Clarity',
      AUC = auc_clarity,
      Cluster = c
    )
    results <- results %>% add_row(
      Seed = t, 
      Dimension = 10, 
      Method = 'Combined',
      AUC = auc_comb,
      Cluster = c
    )    
    
  }
  return(results)
}

write.csv(aucs_cluster, 'simulations/results/cluster_multiway_aucs.csv')
stopCluster(cluster)
registerDoSEQ()
