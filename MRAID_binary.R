# install.packages("invgamma")

library(truncnorm)
library(Rcpp)
library(invgamma)


sourceCpp(code='
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
double safe_log(double x) {
  double epsilon = 1e-8;
  return log(std::max(x, epsilon));;
}
// [[Rcpp::export]]
double safe_exp(double x) {
  if (x > 709) return exp(709);  
  if (x < -745) return exp(-745); 
  return exp(x);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd removeColumn(Eigen::MatrixXd matrix, int j) {
    int nrow = matrix.rows(); 
    int ncol = matrix.cols(); 
    if (j < 0 || j >= ncol) {
        Rcpp::stop("Index j is out of range!");
    }
    Eigen::MatrixXd result(nrow, ncol - 1);
    if (j > 0) {
        result.block(0, 0, nrow, j) = matrix.block(0, 0, nrow, j);
    }
    if (j < ncol - 1) {
        result.block(0, j, nrow, ncol - j - 1) = matrix.block(0, j + 1, nrow, ncol - j - 1);
    }
    return result;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd removeElement(Eigen::VectorXd vec, int j) {
    int n = vec.size(); 
    if (j < 0 || j >= n) {
        Rcpp::stop("Index j is out of range!");
    }
    Eigen::VectorXd result(n - 1);
    for (int i = 0, newIndex = 0; i < n; ++i) {
        if (i == j) continue; 
        result[newIndex++] = vec[i];
    }
    return result;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List update_parameters(Eigen::VectorXd beta_iter, Eigen::VectorXd ita_u_iter, 
                       Eigen::VectorXd gamma_iter, Eigen::VectorXd tau_iter, Eigen::VectorXd v_iter,
                       double alpha_iter, double rho_iter, double sigmaY_iter, double sigmaX_iter, 
                       double sigma_beta_iter, double sigma_ita_iter,
                       Eigen::MatrixXd Zy, Eigen::MatrixXd Zx, Eigen::VectorXd ystar, Eigen::VectorXd x,
                       int n1, int n2, double pi_beta_iter, double pi_1_iter, double pi_0_iter,
                       double pi_c_iter) {
  int J = beta_iter.size(); // Number of variables

  for (int j = 0; j < J; ++j) {
    double K1, mu_beta_j, K2, mu_ita_j;
    double log_prob_gamma_1, log_prob_gamma_0, log_prob_tau_1, log_prob_tau_0, log_prob_v_1, log_prob_v_0, max_log_prob, prob_sum;
    double prob_gamma_1, prob_gamma_0, prob_tau_1, prob_tau_0, prob_v_1, prob_v_0;

    // Update beta_j
    if (gamma_iter(j) == 1) {
      K1 = pow(alpha_iter + rho_iter * v_iter(j), 2) * (n2 - 1) / pow(sigmaY_iter, 2)
           + (n1 - 1) / pow(sigmaX_iter, 2) + 1 / pow(sigma_beta_iter, 2);
      if (K1 < 1e-8) K1 = 1e-8;
      
      Eigen::MatrixXd Zy_sub = removeColumn(Zy, j);
      Eigen::VectorXd beta_sub = removeElement(beta_iter, j);
      Eigen::VectorXd v_sub = removeElement(v_iter, j);
      Eigen::VectorXd residual = ystar - alpha_iter * Zy_sub * beta_sub -
                           rho_iter * Zy_sub * beta_sub.cwiseProduct(v_sub) -
                           Zy * ita_u_iter;
      Eigen::MatrixXd Zx_sub = removeColumn(Zx, j);
      Eigen::VectorXd x_residual = x - Zx_sub * beta_sub;
      double term1 = (residual.transpose() * Zy.col(j)).value() / pow(sigmaY_iter, 2);
      double term2 = (x_residual.transpose() * Zx.col(j)).value() / pow(sigmaX_iter, 2);
      mu_beta_j = ((alpha_iter + rho_iter * v_iter(j)) * term1 + term2) / K1;

      beta_iter(j) = R::rnorm(mu_beta_j, sqrt(1 / K1));
    } else {
      K1 = pow(alpha_iter + rho_iter * v_iter(j), 2) * (n2 - 1) / pow(sigmaY_iter, 2)
           + (n1 - 1) / pow(sigmaX_iter, 2);
      if (K1 < 1e-8) K1 = 1e-8;
      
      Eigen::MatrixXd Zy_sub = removeColumn(Zy, j);
      Eigen::VectorXd beta_sub = removeElement(beta_iter, j);
      Eigen::VectorXd v_sub = removeElement(v_iter, j);
      Eigen::VectorXd residual = ystar - alpha_iter * Zy_sub * beta_sub -
                           rho_iter * Zy_sub * beta_sub.cwiseProduct(v_sub) -
                           Zy * ita_u_iter;
      Eigen::MatrixXd Zx_sub = removeColumn(Zx, j);
      Eigen::VectorXd x_residual = x - Zx_sub * beta_sub;
      double term1 = (residual.transpose() * Zy.col(j)).value() / pow(sigmaY_iter, 2);
      double term2 = (x_residual.transpose() * Zx.col(j)).value() / pow(sigmaX_iter, 2);
      mu_beta_j = ((alpha_iter + rho_iter * v_iter(j)) * term1 + term2) / K1;
      
      beta_iter(j) = 0;
    }
    
    
    // Update ita_uj 
    if (tau_iter(j) == 1) {
      K2 = (n2 - 1) / pow(sigmaY_iter, 2) + 1 / pow(sigma_ita_iter, 2);
      if (K2 < 1e-8) K2 = 1e-8;
      Eigen::VectorXd ita_u_sub = removeElement(ita_u_iter, j);
      Eigen::MatrixXd Zy_sub = removeColumn(Zy, j);
      Eigen::VectorXd residual = ystar - Zy_sub * ita_u_sub -
                           alpha_iter * Zy * beta_iter - 
                           rho_iter*Zy*beta_iter.cwiseProduct(v_iter);
      double term1 = (residual.transpose() * Zy.col(j)).value() / pow(sigmaY_iter, 2);
      mu_ita_j = term1/ K2;
      ita_u_iter(j) = R::rnorm(mu_ita_j, sqrt(1 / K2));
    } else {
      K2 = (n2 - 1) / pow(sigmaY_iter, 2);
      if (K2 < 1e-8) K2 = 1e-8;
      Eigen::VectorXd ita_u_sub = removeElement(ita_u_iter, j);
      Eigen::MatrixXd Zy_sub = removeColumn(Zy, j);
      Eigen::VectorXd residual = ystar - Zy_sub * ita_u_sub -
                           alpha_iter * Zy * beta_iter - 
                           rho_iter*Zy*beta_iter.cwiseProduct(v_iter);
      double term1 = (residual.transpose() * Zy.col(j)).value() / pow(sigmaY_iter, 2);
      mu_ita_j = term1/ K2;
      ita_u_iter(j) = 0;
    }
    
    
    // Update gamma_j
    if (gamma_iter(j) == 1) {
      log_prob_gamma_1 = (pow(mu_beta_j, 2) / (2 * 1 / K1)) + 0.5 * safe_log(1 / K1) - 0.5 * safe_log(pow(sigma_beta_iter, 2)) + safe_log(pi_beta_iter) + tau_iter(j) * safe_log(pi_1_iter) + (1 - tau_iter(j)) * safe_log(1 - pi_1_iter) + v_iter(j) * safe_log(pi_c_iter) + (1 - v_iter(j)) * safe_log(1 - pi_c_iter);
      log_prob_gamma_0 = safe_log(1 - pi_beta_iter) + tau_iter(j) * safe_log(pi_0_iter) + 
                          (1 - tau_iter(j)) * safe_log(1 - pi_0_iter);
  
      double max_log_prob = std::max(log_prob_gamma_1, log_prob_gamma_0);
      prob_gamma_1 = safe_exp(log_prob_gamma_1 - max_log_prob);
      prob_gamma_0 = safe_exp(log_prob_gamma_0 - max_log_prob);
      double prob_sum = prob_gamma_1 + prob_gamma_0;
      gamma_iter(j) = R::rbinom(1, prob_gamma_1 / prob_sum);
    } else {
      log_prob_gamma_1 = (pow(mu_beta_j, 2) / (2 * 1 / K1)) + 0.5 * safe_log(1 / K1) - 0.5 * safe_log(pow(sigma_beta_iter, 2)) + safe_log(pi_beta_iter) + tau_iter(j) * safe_log(pi_1_iter) + (1 - tau_iter(j)) * safe_log(1 - pi_1_iter) + v_iter(j) * safe_log(pi_c_iter) + (1 - v_iter(j)) * safe_log(1 - pi_c_iter);
      log_prob_gamma_0 = safe_log(1 - pi_beta_iter) + tau_iter(j) * safe_log(pi_0_iter) + 
                          (1 - tau_iter(j)) * safe_log(1 - pi_0_iter);
  
      double max_log_prob = std::max(log_prob_gamma_1, log_prob_gamma_0);
      prob_gamma_1 = safe_exp(log_prob_gamma_1 - max_log_prob);
      prob_gamma_0 = safe_exp(log_prob_gamma_0 - max_log_prob);
      double prob_sum = prob_gamma_1 + prob_gamma_0;
      gamma_iter(j) = R::rbinom(1, prob_gamma_1 / prob_sum);
    }
    
    
    // Update tau_j
    if (tau_iter(j) == 1) {
      log_prob_tau_1 = (pow(mu_ita_j, 2) / (2 * 1 / K2)) + 0.5 * safe_log(1 / K2) - 
                        0.5 * safe_log(pow(sigma_ita_iter, 2)) + gamma_iter(j) * safe_log(pi_1_iter) + 
                        (1 - gamma_iter(j)) * safe_log(pi_0_iter);
      log_prob_tau_0 = gamma_iter(j) * safe_log(1 - pi_1_iter) + (1 - gamma_iter(j)) * safe_log(1 - pi_0_iter);
  
      double max_log_prob = std::max(log_prob_tau_1, log_prob_tau_0);
      double prob_tau_1 = safe_exp(log_prob_tau_1 - max_log_prob);
      double prob_tau_0 = safe_exp(log_prob_tau_0 - max_log_prob);
      double prob_sum = prob_tau_1 + prob_tau_0;
      tau_iter(j) = R::rbinom(1, prob_tau_1 / prob_sum);
    } else {
      log_prob_tau_1 = (pow(mu_ita_j, 2) / (2 * 1 / K2)) + 0.5 * safe_log(1 / K2) - 
                        0.5 * safe_log(pow(sigma_ita_iter, 2)) + gamma_iter(j) * safe_log(pi_1_iter) + 
                        (1 - gamma_iter(j)) * safe_log(pi_0_iter);
      log_prob_tau_0 = gamma_iter(j) * safe_log(1 - pi_1_iter) + (1 - gamma_iter(j)) * safe_log(1 - pi_0_iter);
  
      max_log_prob = std::max(log_prob_tau_1, log_prob_tau_0);
      prob_tau_1 = safe_exp(log_prob_tau_1 - max_log_prob);
      prob_tau_0 = safe_exp(log_prob_tau_0 - max_log_prob);
      prob_sum = prob_tau_1 + prob_tau_0;
      tau_iter(j) = R::rbinom(1, prob_tau_1 / prob_sum);
    }
    
    
     // Update v_j
    Eigen::MatrixXd Zy_sub = removeColumn(Zy, j);
    Eigen::VectorXd beta_sub = removeElement(beta_iter, j);
    Eigen::VectorXd v_sub = removeElement(v_iter, j);
    Eigen::VectorXd residual = ystar - alpha_iter * Zy * beta_iter -
                         rho_iter * Zy_sub * beta_sub.cwiseProduct(v_sub) -
                         Zy * ita_u_iter;
    double term1 = (residual.transpose() * Zy.col(j)).value()*beta_iter(j)*rho_iter;
    log_prob_v_1 = -(pow((rho_iter * beta_iter(j)), 2) * (n2 - 1) - 
                      2*term1) / (2 * pow(sigmaY_iter, 2)) + gamma_iter(j) * safe_log(pi_c_iter);
    log_prob_v_0 = gamma_iter(j) * safe_log(1 - pi_c_iter);

    max_log_prob = std::max(log_prob_v_1, log_prob_v_0);
    prob_v_1 = safe_exp(log_prob_v_1 - max_log_prob);
    prob_v_0 = safe_exp(log_prob_v_0 - max_log_prob);
    prob_sum = prob_v_1 + prob_v_0;
    v_iter(j) = R::rbinom(1, prob_v_1 / prob_sum);
  
    
  }    
  return List::create(Named("beta_iter") = beta_iter,
                      Named("ita_u_iter") = ita_u_iter,
                      Named("gamma_iter") = gamma_iter,
                      Named("tau_iter") = tau_iter,
                      Named("v_iter") = v_iter);
}
')


gibbs_probit_binary= function(x, y,Zx,Zy, lambda_beta1, lambda_beta2, lambda_c1, 
                              lambda_c2, lambda_21,lambda_22, lambda_31, 
                              lambda_32,a_beta,b_beta,a_ita,b_ita,
                              n_iter=1000,binaryX=FALSE,binaryY=TRUE,
                              alpha_initial,burnin_proportion=0.2) {
  a = Sys.time()
  n1 = length(x)
  n2 = length(y)
  p = dim(Zx)[2]
  ystar = rep(0, n2)
  xstar = rep(0,n1)
  beta_res = c()
  alpha_res = c()
  ystar_res = c()
  rho_res = c()
  sigma_beta_res = c()
  sigma_ita_res = c()
  sigmaY_res = c()
  sigmaX_res = c()
  
  beta_iter = c(rep(1,round(p/2)),rep(0,p-round(p/2))) 
  gamma_iter = c(rep(1,round(p/2)),rep(0,p-round(p/2))) 
  v_iter = c(rep(1,round(p/2)),rep(0,p-round(p/2))) 
  ita_u_iter = c(rep(1,round(p/2)),rep(0,p-round(p/2))) 
  tau_iter = c(rep(1,round(p/2)),rep(0,p-round(p/2))) 
  alpha_iter = alpha_initial
  rho_iter = 0.5
  sigmaY_iter = sqrt(var(y))
  sigmaX_iter = sqrt(var(x))
  sigma_beta_iter = 1
  sigma_ita_iter = 1
  pi_beta_iter = 0.1
  pi_c_iter = 0.05
  pi_1_iter = 0.25  
  pi_0_iter = 0.005
  for (i in 1:n_iter) {
    ### update the latent variable
    if(binaryY==T){
      mu_ystar =  Zy%*%beta_iter*alpha_iter + Zy%*%(beta_iter*v_iter)*rho_iter +
        Zy%*%ita_u_iter
      ystar[y == 0] =  rtruncnorm(sum(y == 0), a = -Inf, b = 0, 
                                  mean = mu_ystar[y == 0],sd = sigmaY_iter)
      ystar[y == 1] = rtruncnorm(sum(y == 1), a = 0, b = Inf, 
                                 mean = mu_ystar[y == 1],sd = sigmaY_iter)
    }else{
      ystar = y
    }
    
    if(binaryX==T){
      mu_xstar =  Zx%*%beta_iter
      xstar[x == 0] =  rtruncnorm(sum(x == 0), a = -Inf, b = 0, 
                                  mean = mu_xstar[x == 0],sd = sigmaX_iter)
      xstar[x == 1] = rtruncnorm(sum(x == 1), a = 0, b = Inf, 
                                 mean = mu_xstar[x == 1],sd = sigmaX_iter)
    }else{
      xstar = x
    }
    
    # for(j in 1:length(beta_iter)){
    #   ### sample beta given gamma_j
    #   if(gamma_iter[j]==1){ # sample beta given gamma_j = 1
    #     K1 = (alpha_iter+rho_iter*v_iter[j])^2*(n2-1)/sigmaY_iter^2 +
    #       (n1-1)/sigmaX_iter^2+1/sigma_beta_iter^2
    #   # since Zx and Zy are standardized, we have t(Zyj)%*%Zyj=n2-1, t(Zxj)%*%Zxj=n1-1
    #     mu_beta_j = ((alpha_iter+rho_iter*v_iter[j])*
    #                    t(ystar-alpha_iter*Zy[,-j]%*%beta_iter[-j]-
    #                        rho_iter*Zy[,-j]%*%(beta_iter[-j]*v_iter[-j])-
    #                        Zy%*%ita_u_iter)%*%Zy[,j]/sigmaY_iter^2 +
    #                    t(xstar-Zx[,-j]%*%beta_iter[-j])%*%Zx[,j]/sigmaX_iter^2)/K1
    #     beta_iter[j] = rnorm(1,mean=mu_beta_j,sd=sqrt(1/K1))
    #   }else{ # sample beta given gamma_j = 0
    #     K1 = (alpha_iter+rho_iter*v_iter[j])^2*(n2-1)/sigmaY_iter^2 +
    #       (n1-1)/sigmaX_iter^2
    #   # since Zx and Zy are standardized, we have t(Zyj)%*%Zyj=n2-1, t(Zxj)%*%Zxj=n1-1
    #     mu_beta_j = ((alpha_iter+rho_iter*v_iter[j])*
    #                    t(ystar-alpha_iter*Zy[,-j]%*%beta_iter[-j]-
    #                        rho_iter*Zy[,-j]%*%(beta_iter[-j]*v_iter[-j])-
    #                        Zy%*%ita_u_iter)%*%Zy[,j]/sigmaY_iter^2 +
    #                    t(xstar-Zx[,-j]%*%beta_iter[-j])%*%Zx[,j]/sigmaX_iter^2)/K1
    #     beta_iter[j] = 0
    #   }
    # 
    #   ### sample ita_uj given tau_j
    #   if(tau_iter[j]==1){ # sample ita_uj given tau_j = 1
    #     K2 = (n2-1)/sigmaY_iter^2 + 1/sigma_ita_iter^2
    #     # since Zx and Zy are standardized, we have t(Zyj)%*%Zyj=n2-1, t(Zxj)%*%Zxj=n1-1
    #     mu_ita_j = (t(ystar-Zy[,-j]%*%ita_u_iter[-j]-alpha_iter*Zy%*%beta_iter-
    #                     rho_iter*Zy%*%(beta_iter*v_iter))%*%Zy[,j]/sigmaY_iter^2)/K2
    #     ita_u_iter[j] = rnorm(1,mean=mu_ita_j,sd=sqrt(1/K2))
    #   }else{ # sample ita_uj given tau_j = 0
    #     K2 = (n2-1)/sigmaY_iter^2
    #     # since Zx and Zy are standardized, we have t(Zyj)%*%Zyj=n2-1, t(Zxj)%*%Zxj=n1-1
    #     mu_ita_j = (t(ystar-Zy[,-j]%*%ita_u_iter[-j]-alpha_iter*Zy%*%beta_iter-
    #                     rho_iter*Zy%*%(beta_iter*v_iter))%*%Zy[,j]/sigmaY_iter^2)/K2
    #     ita_u_iter[j] = 0
    #   }
    # 
    #   ### sample gamma_j
    #   log_prob_gamma_1 = (mu_beta_j^2 / (2*1/K1)) + 0.5*log(1/K1) -
    #     0.5*log(sigma_beta_iter^2) + log(pi_beta_iter)+
    #     tau_iter[j] * log(pi_1_iter) + (1 - tau_iter[j]) * log(1 - pi_1_iter)+
    #     v_iter[j] * log(pi_c_iter)+ (1 - v_iter[j]) * log(1 - pi_c_iter)
    #   log_prob_gamma_0 <- log(1-pi_beta_iter) + tau_iter[j] * log(pi_0_iter) +
    #                       (1 - tau_iter[j]) * log(1 - pi_0_iter)
    #   if (is.infinite(exp(log_prob_gamma_1)) || is.infinite(exp(log_prob_gamma_0))) {
    #     cat("Detected Inf in gamma_j probabilities at j =", j,
    #         ", applying Log-Sum-Exp technique...\n")
    #     # stablization
    #     max_log_prob <- max(log_prob_gamma_1, log_prob_gamma_0)
    #     prob_gamma_1 <- exp(log_prob_gamma_1 - max_log_prob)
    #     prob_gamma_0 <- exp(log_prob_gamma_0 - max_log_prob)
    #   }else{
    #     prob_gamma_1 <- exp(log_prob_gamma_1)
    #     prob_gamma_0 <- exp(log_prob_gamma_0)
    #   }
    #   prob_sum <- prob_gamma_1 + prob_gamma_0
    #   prob_gamma_1 <- prob_gamma_1 / prob_sum # normalizated prob
    #   prob_gamma_0 <- prob_gamma_0 / prob_sum # normalizated prob
    #   gamma_iter[j] <- rbinom(1, size = 1, prob = prob_gamma_1)
    # 
    # 
    #   ### sample tau_j
    #   log_prob_tau_1 = (mu_ita_j^2 / (2*1/K2)) + 0.5*log(1/K2) -
    #     0.5*log(sigma_ita_iter^2) + gamma_iter[j]*log(pi_1_iter)+
    #     (1-gamma_iter[j])*log(pi_0_iter)
    #   log_prob_tau_0 <- gamma_iter[j]*log(1-pi_1_iter) +
    #                       (1 - gamma_iter[j])*log(1 - pi_0_iter)
    #   if (is.infinite(exp(log_prob_tau_1)) || is.infinite(exp(log_prob_tau_0))) {
    #     cat("Detected Inf in tau_j probabilities at j =", j,
    #         ", applying Log-Sum-Exp technique...\n")
    #     # stablization
    #     max_log_prob <- max(log_prob_tau_1, log_prob_tau_0) # 提取最大 log 值
    #     prob_tau_1 <- exp(log_prob_tau_1 - max_log_prob)
    #     prob_tau_0 <- exp(log_prob_tau_0 - max_log_prob)
    #   }else{
    #     prob_tau_1 <- exp(log_prob_tau_1)
    #     prob_tau_0 <- exp(log_prob_tau_0)
    #   }
    #   prob_sum <- prob_tau_1 + prob_tau_0
    #   prob_tau_1 <- prob_tau_1 / prob_sum # normalizated prob
    #   prob_tau_0 <- prob_tau_0 / prob_sum # normalizated prob
    #   tau_iter[j] <- rbinom(1, size = 1, prob = prob_tau_1)
    # 
    #   ### sample v_j
    #   log_prob_v_1 = -((rho_iter^2*beta_iter[j]^2*(n2-1)-
    #                      2*t(ystar-alpha_iter*Zy%*%beta_iter-
    #                           rho_iter*Zy[,-j]%*%(beta_iter[-j]*v_iter[-j])-
    #                           Zy%*%ita_u_iter)%*%Zy[,j]*beta_iter[j]*rho_iter)/
    #                      (2*sigmaY_iter^2))+gamma_iter[j]*log(pi_c_iter)
    #   # since Zx and Zy are standardized, we have t(Zyj)%*%Zyj=n2-1, t(Zxj)%*%Zxj=n1-1
    #   log_prob_v_0 <- gamma_iter[j]*log(1-pi_c_iter)
    #   if (is.infinite(exp(log_prob_v_1)) || is.infinite(exp(log_prob_v_0))) {
    #     cat("Detected Inf in v_j probabilities at j =", j,
    #         ", applying Log-Sum-Exp technique...\n")
    #     # stablization
    #     max_log_prob <- max(log_prob_v_1, log_prob_v_0)
    #     prob_v_1 <- exp(log_prob_v_1 - max_log_prob)
    #     prob_v_0 <- exp(log_prob_v_0 - max_log_prob)
    #   }else{
    #     prob_v_1 <- exp(log_prob_v_1)
    #     prob_v_0 <- exp(log_prob_v_0)
    #   }
    #   prob_sum <- prob_v_1 + prob_v_0
    #   prob_v_1 <- prob_v_1 / prob_sum # normalizated prob
    #   prob_v_0 <- prob_v_0 / prob_sum # normalizated prob
    #   v_iter[j] <- rbinom(1, size = 1, prob = prob_v_1)
    # }
    
    ### update posterior j-dim parameters
    result <- update_parameters(beta_iter = beta_iter, ita_u_iter=ita_u_iter,
                                gamma_iter=gamma_iter, tau_iter=tau_iter, v_iter=v_iter,
                                alpha_iter=alpha_iter,rho_iter=rho_iter,
                                sigmaY_iter=sigmaY_iter, sigmaX_iter=sigmaX_iter,
                                sigma_beta_iter=sigma_beta_iter,
                                sigma_ita_iter=sigma_ita_iter,
                                Zy=Zy, Zx=Zx, ystar=ystar, x=xstar, n1=n1, n2=n2,
                                pi_beta_iter=pi_beta_iter,
                                pi_1_iter=pi_1_iter, pi_0_iter=pi_0_iter,
                                pi_c_iter=pi_c_iter)
    beta_iter = result$beta_iter
    ita_u_iter = result$ita_u_iter
    gamma_iter = result$gamma_iter
    tau_iter = result$tau_iter
    v_iter = result$v_iter
    
    ### update posterior rho
    var_rho = (t(v_iter*beta_iter)%*%t(Zy)%*%Zy%*%(v_iter*beta_iter)/sigmaY_iter^2)^(-1)
    mean_rho = (t(v_iter*beta_iter)%*%t(Zy)%*%ystar-t(v_iter*beta_iter)%*%t(Zy)%*%Zy%*%
                  ita_u_iter-t(v_iter*beta_iter)%*%t(Zy)%*%Zy%*%(beta_iter)*alpha_iter
    )/sigmaY_iter^2*var_rho
    rho_iter = rnorm(1,mean=mean_rho,sd=sqrt(var_rho))
    rho_iter <- pmax(-1, pmin(1, rho_iter)) ##########set rho btw -1 and 1
    
    ### update posterior alpha
    var_alpha = (t(beta_iter)%*%t(Zy)%*%Zy%*%beta_iter/sigmaY_iter^2)^(-1)
    mean_alpha = (t(beta_iter)%*%t(Zy)%*%ystar-t(beta_iter)%*%t(Zy)%*%Zy%*%
                    ita_u_iter-t(beta_iter)%*%t(Zy)%*%Zy%*%(beta_iter*v_iter)*rho_iter
    )/sigmaY_iter^2*var_alpha
    alpha_iter = rnorm(1,mean=mean_alpha,sd=sqrt(var_alpha))
    
    
    ### update posterior pi_beta, pi_1, pi_0, pi_c
    pi_beta_iter = rbeta(1, shape1 = sum(gamma_iter)+lambda_beta1, 
                         shape2 =  sum(1-gamma_iter)+lambda_beta2)
    pi_1_iter = rbeta(1, shape1 = sum(gamma_iter*tau_iter)+lambda_21, 
                      shape2 =  sum(gamma_iter*(1-tau_iter))+lambda_22)
    pi_0_iter = rbeta(1, shape1 = sum((1-gamma_iter)*tau_iter)+lambda_31, 
                      shape2 =  sum((1-gamma_iter)*(1-tau_iter))+lambda_32)
    pi_c_iter = rbeta(1, shape1 = sum(gamma_iter*v_iter)+lambda_c1, 
                      shape2 =  sum(gamma_iter*(1-v_iter))+lambda_c1)
    
    ### update posterior sigmaX
    shape_param_sigmaX = n1/2-1
    scale_param_sigmaX = (t(xstar)%*%xstar+t(beta_iter)%*%t(Zx)%*%Zx%*%beta_iter-
                            2*t(beta_iter)%*%t(Zx)%*%xstar)/2
    sigmaX_iter = sqrt(rinvgamma(1, shape = shape_param_sigmaX, 
                                 scale = as.numeric(scale_param_sigmaX)))
    ### update posterior sigmaY
    shape_param_sigmaY = n2/2-1
    scale_param_sigmaY = t(ystar - alpha_iter*Zy%*%beta_iter -
                             rho_iter*Zy%*%(beta_iter*v_iter)-Zy%*%ita_u_iter) %*%
      (ystar - alpha_iter*Zy%*%beta_iter -rho_iter*Zy%*%(beta_iter*v_iter)-
         Zy%*%ita_u_iter)/2
    sigmaY_iter = sqrt(rinvgamma(1, shape = shape_param_sigmaY, 
                                 scale = as.numeric(scale_param_sigmaY)))
    ### update posterior sigma_beta
    shape_param_sigma_beta = sum(gamma_iter)/2+a_beta
    scale_param_sigma_beta = sum(gamma_iter*beta_iter^2)/2+b_beta
    sigma_beta_iter = sqrt(rinvgamma(1, shape = shape_param_sigma_beta, 
                                     scale = as.numeric(scale_param_sigma_beta)))
    ### update posterior sigma_ita
    shape_param_sigma_ita = sum(tau_iter)/2+a_ita
    scale_param_sigma_ita = sum(tau_iter*ita_u_iter^2)/2+b_ita
    sigma_beta_iter = sqrt(rinvgamma(1, shape = shape_param_sigma_ita, 
                                     scale = as.numeric(scale_param_sigma_ita)))
    
    beta_res = rbind(beta_res, beta_iter)
    ystar_res = rbind(ystar_res, ystar)
    alpha_res = c(alpha_res, alpha_iter)
    rho_res = c(rho_res,rho_iter)
    sigma_beta_res = c(sigma_beta_res,sigma_beta_iter)
    sigma_ita_res = c(sigma_ita_res,sigma_ita_iter)
    sigmaX_res = c(sigmaX_res,sigmaX_iter)
    sigmaY_res = c(sigmaY_res,sigmaY_iter)
    # print(i)
  }
  b = Sys.time()
  print(b - a)
  burnin = round(burnin_proportion*n_iter)
  
  beta_est = mean(beta_res[burnin:n_iter])
  #ystar_est = colmeans(ystar_res[burnin:n_iter,])
  alpha_est = mean(alpha_res[burnin:n_iter])
  alpha_sd = sd(alpha_res[burnin:n_iter])
  rho_est = mean(rho_res[burnin:n_iter])
  sigma_beta_est = mean(sigma_beta_res[burnin:n_iter])
  sigma_ita_est = mean(sigma_ita_res[burnin:n_iter])
  sigmaX_est = mean(sigmaX_res[burnin:n_iter])
  sigmaY_est = mean(sigmaY_res[burnin:n_iter])
  return(list(alpha = alpha_est, rho = rho_est, sigma_beta = sigma_beta_est, 
              alpha_sd = alpha_sd,sigma_ita = sigma_ita_est, sigmaX = sigmaX_est,
              sigmaY = sigmaY_est,beta_res = beta_res, ystar_res = ystar_res, 
              alpha_res = alpha_res))
}


# Function to rerun the target function immediately when specific warnings occur
run_with_warning_check <- function(retry_function, max_attempts = 15, ..., warning_keyword = "NAs produced") {
  # retry_function: The target function to execute (e.g., gibbs_probit_binary)
  # max_attempts: Maximum number of retry attempts
  # warning_keyword: The specific warning message to trigger a retry
  # ...: Arguments to pass to the target function
  
  attempt <- 1  # Initialize attempt counter
  result <- NULL  # Variable to store function result
  
  repeat {
    # Try running the function and handle warnings
    result <- tryCatch(
      withCallingHandlers(
        retry_function(...),  # Execute the target function
        warning = function(w) {
          # Check if the warning message contains the specified keyword
          if (grepl(warning_keyword, conditionMessage(w))) {
            message(sprintf("Warning detected on attempt %d: %s", attempt, conditionMessage(w)))
            stop("Triggering function retry due to specific warning.")  # Terminate the current execution
          }
        }
      ),
      error = function(e) {
        # Handle errors (e.g., triggered by stop)
        message(sprintf("Retrying function due to warning on attempt %d...", attempt))
        NULL
      }
    )
    
    # Exit if successful or maximum attempts reached
    if (!is.null(result) || attempt >= max_attempts) {
      if (is.null(result) && attempt >= max_attempts) {
        message("Maximum attempts reached. Returning NULL.")
      } else {
        message("Function completed successfully without critical warnings.")
      }
      break
    }
    
    # Increment attempt counter and retry
    attempt <- attempt + 1
  }
  
  return(result)  # Return the result or NULL if failed
}
