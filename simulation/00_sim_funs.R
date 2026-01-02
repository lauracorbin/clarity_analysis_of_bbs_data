sim_data <- function(
    n = 500,                      # number of samples observations
    nU = 2,                       # number of unobserved confounders
    nC = 3,                       # number of observed confounders
    nX = 100,                     # number of X variables
    clusters = 6,                 # number of clusters
    rho_time = 0.7,               # correlation of same coefficient across times
    sigma_x = 0.3,                # covariance for effects outside of activated clusters
    sigma_c = 0.7,                # covariance for effects within activated clusters
    mean_effect = 1,              # mean of pathway C*
    covar_means = NULL,           # numeric or vector mean of covariates 
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  
  # ---------- 1. Generate correlated confounders (U + C) ----------
  pUC <- nU + nC
  if (is.null(covar_means)) {
    means <- runif(pUC,0,0.5)
  } else if (length(covar_means) == pUC){
    means <- covar_means
  } else if (length(covar_means) == 1) {
    means <- rep(covar_means, pUC)
  } else {
    return(paste("Covar means must of length 1 or nU + nC"))
  }
  
  Sigma_UC <- matrix(0.5, nrow = pUC, ncol = pUC)
  diag(Sigma_UC) <- 1
  UC <- MASS::mvrnorm(n = n, mu = means, Sigma = Sigma_UC)
  colnames(UC) <- c(paste0("U", 1:nU), paste0("C", 1:nC))
  UC_mat <- as.matrix(UC) # n x pUC
  
  # Extra confounder (affects subset of Xs at t2)
  Cstar <- rnorm(n, mean = mean_effect, sd = 0.5) 
  
  
  # ---------- 2. Create clusters structure ----------
  
  groups <- cut(1:nX, clusters, labels = FALSE)
  
  
  # ---------- 3. Draw coefficients for linear model ----------
  
  
  p_block <- nX
  p_total <- 2 * p_block
  
  
  Sigma_base <- matrix(sigma_x, nrow = p_block, ncol= p_block)
  
  coefs_t1 <- matrix(0, nrow = nX, ncol = pUC)
  coefs_t2 <- matrix(0, nrow = nX, ncol = pUC)
  act_clusters <- matrix(0, nrow = (nU + nC), ncol =2)
  
  for (k in 1:pUC) {
    mu_coef <- rep(0, p_total)
    npaths <- 0
    Sigma_C <- Sigma_base 
    paths <- sample(1:clusters, 2)
    act_clusters[k,] <- paths
    inds <- groups %in% paths
    Sigma_C[inds, inds] <- (sigma_x + sigma_c) / 2
    for (p in paths){
      Sigma_C[groups == p, groups == p] <- sigma_c
    }
    
    Sigma_coef <- rbind(
      cbind(Sigma_C, rho_time * Sigma_C),
      cbind(rho_time * Sigma_C, Sigma_C)
    ) 
    
    draw <- MASS::mvrnorm(n = 1, mu = mu_coef, Sigma = Sigma_coef)
    coefs_t1[, k] <- draw[1:p_block]
    coefs_t2[, k] <- draw[(p_block + 1):p_total]
  }
  colnames(coefs_t1) <- colnames(coefs_t2) <- c(paste0("U", 1:nU), paste0("C", 1:nC))
  
  # ---------- 4. Decide which Xs get Cstar at t2 ----------
  
  anomaly_groups <- sample(1:clusters,2)
  
  anomaly_group1 <- anomaly_groups[1]
  anomaly_group2 <- anomaly_groups[2]
  anomalies1_group <- which(groups == anomaly_group1)
  anomalies1 <- anomalies1_group[1:round(length(anomalies1_group)/2)]
  
  #added
  anomalies3 <- anomalies1_group[round(length(anomalies1_group)/2):length(anomalies1_group)]
  
  anomalies2 <- which(groups == anomaly_group2)
  anomalies2 <- anomalies2[1:round(length(anomalies2))]
  
  Cstar_coefs <- rep(0, nX)
  Cstar_coefs[c(anomalies1,anomalies2)] <- runif(length(c(anomalies1, anomalies2)), min = 0.8, max = 1)  
  
  #anomalies <- numeric(nX)
  #anomalies[anomalies1] <- 1
  #anomalies[anomalies2] <- 2
  
  anomalies <- numeric(nX)
  anomalies[anomalies1] <- 2
  anomalies[anomalies2] <- 1
  anomalies[anomalies3] <- 3
  
  # ---------- 5. Create X errors ----------
  
  Sigma_X <- matrix(0, nX, nX)
  diag(Sigma_X) <- 1
  
  E1 <- MASS::mvrnorm(n = n, mu = rep(0, nX), Sigma = Sigma_X)
  E2 <- MASS::mvrnorm(n = n, mu = rep(0, nX), Sigma = Sigma_X)
  
  # ---------- 6. Build X1 and X2 ----------
  X1 <- matrix(NA, nrow = n, ncol = nX)
  X2 <- matrix(NA, nrow = n, ncol = nX)
  for (k in 1:nX) {
    # linear predictor from U and C at t1
    X1[, k] <- UC_mat %*% as.numeric(coefs_t1[k, ]) + E1[, k] 
    # t2: AR term + UC effects (t2 coefs) + Cstar if active + error
    X2[, k] <- #ar_coef * X1[, k] +
      UC_mat %*% as.numeric(coefs_t2[k, ]) +
      Cstar * Cstar_coefs[k] + 
      E2[, k]
  }
  colnames(X1) <- paste0("X", 1:nX)
  colnames(X2) <- paste0("X", 1:nX)
  
  # ---------- 7. Pack effect-size tables ----------
  # tidy data frames listing coefs for each X
  effects_t1 <- data.frame(X = paste0("X", 1:nX), coefs_t1)
  effects_t2 <- data.frame(X = paste0("X", 1:nX), coefs_t2)
  effects_t2$Cstar <- Cstar_coefs
  
  
  list(
    X1 = X1,
    X2 = X2,
    UC = UC_mat,
    effects_t1 = effects_t1,
    effects_t2 = effects_t2,
    anomalies = anomalies,
    act_clusters = act_clusters
  )
  
} 

  

