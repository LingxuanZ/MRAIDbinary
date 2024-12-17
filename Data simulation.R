library(Rcpp)
library(MASS)
library(Rfast)

simulate_exposure_outcome <- function(
    n1 = 3000,                 # Number of individuals for exposure data
    n2 = 4000,                 # Number of individuals for outcome data
    p = 6000,                  # Total number of SNPs
    K_number = 100,            # Number of causal SNPs
    PVE_Zx = 0.05,             # Proportion of variance in exposure explained by SNPs
    PVE_alpha = 0.05,          # Proportion of variance in outcome explained by causal effect
    PVE_u = 0.025,             # Proportion of variance in outcome explained by uncorrelated pleiotropy
    rho_correlated_pleiotropy = sqrt(0.02), # Correlation factor for correlated pleiotropy
    pi_c = 0.3,                # Proportion of causal SNPs exhibiting correlated pleiotropy
    pi_1 = 0.1,                # Proportion of causal SNPs exhibiting uncorrelated pleiotropy
    total_uncorrelated_pleiotropy = 100     # Total SNPs with uncorrelated pleiotropy
) {
  # Step 1: Simulate Zx (genotypes for exposure)
  mean_vector <- rep(0, p)
  cov_matrix <- diag(1, p)
  Zx <- Rfast::rmvnorm(n1, mean_vector, sigma = cov_matrix)
  
  # Causal selected SNPs indicator
  K_indicator <- rep(0, p)
  K_indicator[sample(p, K_number)] <- 1
  
  # Extract causal SNPs and calculate beta
  Zx_tilde <- Zx[, which(K_indicator == 1)]
  sd_beta <- sqrt(PVE_Zx / K_number)
  beta <- rep(0, p)
  beta_causal <- rnorm(K_number, mean = 0, sd = sd_beta)
  beta[which(K_indicator == 1)] <- beta_causal
  
  # Simulate exposure (X)
  epsilon_x <- rnorm(n1, mean = 0, sd = sqrt(1 - PVE_Zx))
  X <- Zx_tilde %*% beta_causal + epsilon_x
  X_b <- ifelse(X >= 0, 1, 0)
  
  # Step 2: Simulate Zy (genotypes for outcome)
  Zy <- Rfast::rmvnorm(n2, mean_vector, sigma = cov_matrix)
  Zy_tilde <- Zy[, which(K_indicator == 1)]
  
  # Causal effect
  alpha <- sqrt(PVE_alpha / PVE_Zx)
  
  # Correlated pleiotropy
  correlated_pleiotropy_indicator <- rep(0, p)
  correlated_pleiotropy_indicator[sample(which(K_indicator == 1), round(pi_c * K_number))] <- 1
  Zy_correlated <- Zy[, which(correlated_pleiotropy_indicator == 1)]
  beta_correlated <- beta[which(correlated_pleiotropy_indicator == 1)]
  
  # Uncorrelated pleiotropy
  uncorrelated_causal_count <- round(pi_1 * K_number)
  uncorrelated_noncausal_count <- total_uncorrelated_pleiotropy - uncorrelated_causal_count
  uncorrelated_pleiotropy_indicator <- rep(0, p)
  uncorrelated_pleiotropy_indicator[sample(which(K_indicator == 1), uncorrelated_causal_count)] <- 1
  uncorrelated_pleiotropy_indicator[sample(which(K_indicator == 0), uncorrelated_noncausal_count)] <- 1
  Zy_uncorrelated <- Zy[, which(uncorrelated_pleiotropy_indicator == 1)]
  eta_u <- rnorm(total_uncorrelated_pleiotropy, mean = 0, sd = sqrt(PVE_u / total_uncorrelated_pleiotropy))
  
  # Noise for outcome
  PVE_c <- rho_correlated_pleiotropy^2 * sd_beta^2
  epsilon_y <- rnorm(n2, mean = 0, sd = sqrt(1 - PVE_alpha - PVE_u - PVE_c))
  
  # Simulate outcome (Y)
  Y <- Zy_tilde %*% (alpha * beta_causal) +
    Zy_correlated %*% (rho_correlated_pleiotropy * beta_correlated) +
    Zy_uncorrelated %*% eta_u +
    epsilon_y
  Y_b <- ifelse(Y >= 0, 1, 0)
  
  # Save results
  results <- list(
    X = X,
    X_b = X_b,
    Y = Y,
    Y_b = Y_b,
    Zx = Zx,
    Zy = Zy,
    beta = beta,
    alpha = alpha,
    correlated_pleiotropy_indicator = correlated_pleiotropy_indicator,
    uncorrelated_pleiotropy_indicator = uncorrelated_pleiotropy_indicator
  )
  return(results)
}

results <- simulate_exposure_outcome(n1=300 , n2=400, p=300,PVE_Zx=0.05,PVE_alpha=0.5,rho_correlated_pleiotropy=0,pi_c=0,pi_1=0,total_uncorrelated_pleiotropy=30)
# save(results, file = "simulation1.RData")
list2env(results, envir = .GlobalEnv)
save(list = names(results), file = "simulation_noHorizontal.RData")

results <- simulate_exposure_outcome(n1=300 , n2=400, p=300,PVE_Zx=0.05,PVE_alpha=0.5,rho_correlated_pleiotropy=0,pi_c=0,pi_1=0,total_uncorrelated_pleiotropy=30)
# save(results, file = "simulation1.RData")
list2env(results, envir = .GlobalEnv)
save(list = names(results), file = "simulation_noHorizontal.RData")


results <- simulate_exposure_outcome(n1=300 , n2=400, 
                                     p=300,PVE_Zx=0.05,PVE_u = 0.05,
                                     PVE_alpha=0.5,
                                     rho_correlated_pleiotropy=0,
                                     pi_c=0,pi_1=0.3,total_uncorrelated_pleiotropy=30)
# save(results, file = "simulation1.RData")
list2env(results, envir = .GlobalEnv)
save(list = names(results), file = "simulation_noCorr_wUncorr.RData")

results <- simulate_exposure_outcome(n1=300 , n2=400, 
                                     p=300,PVE_Zx=0.1,PVE_u = 
                                       0.05,PVE_alpha=0.5,
                                     rho_correlated_pleiotropy=sqrt(0.05),
                                     pi_c=0.05,pi_1=0.2,
                                     total_uncorrelated_pleiotropy=30)
# save(results, file = "simulation1.RData")
list2env(results, envir = .GlobalEnv)
save(list = names(results), file = "simulation_wCorrUncorr.RData")
