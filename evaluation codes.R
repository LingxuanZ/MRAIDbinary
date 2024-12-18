####################################################
# no Horizontal pleiotropy
####################################################
PVE_alpha_seq <- seq(0, 0.0025, length.out = 5) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1 = 300, n2 = 400, p = 300, 
                                       PVE_Zx = PVE_Zx, 
                                       PVE_alpha = PVE_alpha,
                                       rho_correlated_pleiotropy = 0, 
                                       pi_c = 0, pi_1 = 0, 
                                       total_uncorrelated_pleiotropy = 30)
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  p = dim(X)[2]
  a_beta = p/10+1
  b_beta = 0.2 # to ensure the prior mean of sigma_beta^2 = 0.2
  a_ita = p/5+1
  b_ita = 0.2 # to ensure the prior mean of sigma_ita^2 = 0.2
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate the 95% confidence intervals for causal effects
ci_upper <- causal_effects + 1.96 * alpha_sds
ci_lower <- causal_effects - 1.96 * alpha_sds

# Determine the y-axis limits
y_min <- min(ci_lower, true_alpha) - 0.01
y_max <- max(ci_upper, true_alpha) + 0.01

# Plot causal effect vs PVE_alpha
par(mar = c(4, 4, 2, 2))  # Adjust margins: bottom, left, top, right
plot(PVE_alpha_seq, causal_effects, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "PVE_alpha", ylab = "Causal Effect Estimates", 
     ylim = c(y_min, y_max),
     main = "Causal Effect vs PVE_alpha with Confidence Intervals")

# Add confidence intervals
lines(PVE_alpha_seq, ci_upper, col = "red", lty = 2)
lines(PVE_alpha_seq, ci_lower, col = "red", lty = 2)

# Add true causal effect line
lines(PVE_alpha_seq, true_alpha, col = "green", lwd = 2, lty = 3)

# Add legend
legend("topright", legend = c("Causal Effect Estimates", "95% CI", "True Causal Effect"), 
       col = c("blue", "red", "green"), lty = c(1, 2, 3), lwd = 2, pch = c(16, NA, NA), 
       bty = "n")


# QQ plot
PVE_alpha_seq <- seq(0, 0.0025, length.out = 30) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1 = 300, n2 = 400, p = 300, 
                                       PVE_Zx = PVE_Zx, 
                                       PVE_alpha = PVE_alpha,
                                       rho_correlated_pleiotropy = 0, 
                                       pi_c = 0, pi_1 = 0, 
                                       total_uncorrelated_pleiotropy = 30)
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate z-scores and p-values
z_scores <- causal_effects / alpha_sds
qqnorm(z_scores)
qqline(z_scores,col="red")


####################################################
# simulations w/o correlated horizontal pleiotropic effect but w/ uncorrelated horizontal pleiotropic effect (PVEu=5%, total # SNPS with uncorrelated pleiotropy=30)
####################################################
PVE_alpha_seq <- seq(0, 0.0025, length.out = 5) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1 = 300, n2 = 400, p = 300, 
                                       PVE_Zx = PVE_Zx, 
                                       PVE_alpha = PVE_alpha,
                                       rho_correlated_pleiotropy = 0, 
                                       pi_c = 0, pi_1 = 0.3, 
                                       total_uncorrelated_pleiotropy = 30)
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate the 95% confidence intervals for causal effects
ci_upper <- causal_effects + 1.96 * alpha_sds
ci_lower <- causal_effects - 1.96 * alpha_sds

# Determine the y-axis limits
y_min <- min(ci_lower, true_alpha) - 0.01
y_max <- max(ci_upper, true_alpha) + 0.01

# Plot causal effect vs PVE_alpha
par(mar = c(4, 4, 2, 2))  # Adjust margins: bottom, left, top, right
plot(PVE_alpha_seq, causal_effects, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "PVE_alpha", ylab = "Causal Effect Estimates", 
     ylim = c(y_min, y_max),
     main = "Causal Effect vs PVE_alpha with Confidence Intervals")

# Add confidence intervals
lines(PVE_alpha_seq, ci_upper, col = "red", lty = 2)
lines(PVE_alpha_seq, ci_lower, col = "red", lty = 2)

# Add true causal effect line
lines(PVE_alpha_seq, true_alpha, col = "green", lwd = 2, lty = 3)

# Add legend
legend("topright", legend = c("Causal Effect Estimates", "95% CI", "True Causal Effect"), 
       col = c("blue", "red", "green"), lty = c(1, 2, 3), lwd = 2, pch = c(16, NA, NA), 
       bty = "n")

# QQ plot
PVE_alpha_seq <- seq(0, 0.0025, length.out = 30) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1 = 300, n2 = 400, p = 300, 
                                       PVE_Zx = PVE_Zx, 
                                       PVE_alpha = PVE_alpha,
                                       rho_correlated_pleiotropy = 0, 
                                       pi_c = 0, pi_1 = 0.3, 
                                       total_uncorrelated_pleiotropy = 30)
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate z-scores and p-values
z_scores <- causal_effects / alpha_sds
qqnorm(z_scores)
qqline(z_scores,col="red")


####################################################
# simulations with both correlated and uncorrelated pleiotropy effect(PVE_u = 0.05,PVE_alpha=0.5, rho_correlated_pleiotropy=sqrt(0.05), pi_c=0.05,pi_1=0.2,total_uncorrelated_pleiotropy=30
####################################################
PVE_alpha_seq <- seq(0, 0.0025, length.out = 5) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1=300 , n2=400, 
                                       p=300,PVE_Zx=0.1,PVE_u = 
                                         0.05,PVE_alpha=0.5,
                                       rho_correlated_pleiotropy=sqrt(0.05),
                                       pi_c=0.05,pi_1=0.2,
                                       total_uncorrelated_pleiotropy=30)
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate the 95% confidence intervals for causal effects
ci_upper <- causal_effects + 1.96 * alpha_sds
ci_lower <- causal_effects - 1.96 * alpha_sds

# Determine the y-axis limits
y_min <- min(ci_lower, true_alpha) - 0.01
y_max <- max(ci_upper, true_alpha) + 0.01

# Plot causal effect vs PVE_alpha
par(mar = c(4, 4, 2, 2))  # Adjust margins: bottom, left, top, right
plot(PVE_alpha_seq, causal_effects, type = "b", pch = 16, col = "blue", lwd = 2,
     xlab = "PVE_alpha", ylab = "Causal Effect Estimates", 
     ylim = c(y_min, y_max),
     main = "Causal Effect vs PVE_alpha with Confidence Intervals")

# Add confidence intervals
lines(PVE_alpha_seq, ci_upper, col = "red", lty = 2)
lines(PVE_alpha_seq, ci_lower, col = "red", lty = 2)

# Add true causal effect line
lines(PVE_alpha_seq, true_alpha, col = "green", lwd = 2, lty = 3)

# Add legend
legend("topright", legend = c("Causal Effect Estimates", "95% CI", "True Causal Effect"), 
       col = c("blue", "red", "green"), lty = c(1, 2, 3), lwd = 2, pch = c(16, NA, NA), 
       bty = "n")

# qq plot
PVE_alpha_seq <- seq(0, 0.0025, length.out = 30) # Set up the PVE_alpha sequence
causal_effects <- numeric(length(PVE_alpha_seq))
alpha_sds <- numeric(length(PVE_alpha_seq))
PVE_Zx = 0.05
true_alpha <- sqrt(PVE_alpha_seq / PVE_Zx)

# Run simulations and store results
for (i in 1:length(PVE_alpha_seq)) {
  PVE_alpha <- PVE_alpha_seq[i]
  
  # Simulate the data
  results <- simulate_exposure_outcome(n1=300 , n2=400, 
                                       p=300,PVE_Zx=0.1,PVE_u = 
                                         0.05,PVE_alpha=0.5,
                                       rho_correlated_pleiotropy=sqrt(0.05),
                                       pi_c=0.05,pi_1=0.2,
                                       total_uncorrelated_pleiotropy=30)
  
  list2env(results, envir = .GlobalEnv)
  
  # Fit initial glm
  n1 <- dim(X)[1]
  fit <- glm(Y_b[1:n1] ~ -1 + X, family = binomial)
  alpha_initial <- abs(summary(fit)$coefficients[1, 1])
  
  # Run Gibbs sampling
  result <- gibbs_probit_binary(x = X, y = Y_b, Zx, Zy, lambda_beta1, lambda_beta2, lambda_c1,
                                binaryX = FALSE, binaryY = TRUE, lambda_c2 = lambda_c2, 
                                lambda_21 = lambda_21, lambda_22 = lambda_22, 
                                lambda_31 = lambda_31, lambda_32 = lambda_32, 
                                a_beta = a_beta, b_beta = b_beta, a_ita = a_ita, b_ita = b_ita, 
                                n_iter = 1000, alpha_initial = alpha_initial, 
                                burnin_proportion = 0.2)
  
  # Store causal effect and its SD
  causal_effects[i] <- result$alpha
  alpha_sds[i] <- result$alpha_sd
}

# Calculate z-scores and p-values
z_scores <- causal_effects / alpha_sds
qqnorm(z_scores)
qqline(z_scores,col="red")
