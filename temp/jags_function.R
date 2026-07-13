
### federated models ==============
FP_continuous <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])
  
  theta[(shared_p + 1)] <- log_sigma2
  
  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha) 
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  # hetero intercept
  for (pp in 1:n_site) {
    theta[(shared_p+1+pp)] <- beta0 + delta[pp]
    delta[pp] ~ dnorm(0, tau_delta[pp])
    tau_delta[pp] <- ifelse(pp == target_site, 1.0E10, 1 / (tau2 * lambda2[pp]))
    lambda2[pp] <- pow(lambda[pp], 2)
    lambda[pp] ~ dt(0, 1, 1) T(0,)
  }
  tau2 <- pow(tau, 2)
  tau ~ dt(0, 1, 1) T(0,)
  
  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.04)

}
"


FP_binary <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha) 
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  # hetero intercept
  for (pp in 1:n_site) {
    theta[(shared_p + pp)] <- beta0 + delta[pp]
    delta[pp] ~ dnorm(0, tau_delta[pp])
    tau_delta[pp] <- ifelse(pp == target_site, 1.0E10, 1 / (tau2 * lambda2[pp]))
    lambda2[pp] <- pow(lambda[pp], 2)
    lambda[pp] ~ dt(0, 1, 1) T(0,)
  }
  tau2 <- pow(tau, 2)
  tau ~ dt(0, 1, 1) T(0,)
  
  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

  ### g-computation in target site
  for (i in 1:N_target) {
    eta[i] <- theta[shared_p + target_site] + inprod(theta[2:beta_p], X_target[i, 2:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


FP_hetero_continuous <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])
  
  theta[(shared_p + 1)] <- log_sigma2
  
  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha) 
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  # site-specific intercepts
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + 1 + pp)] <- beta0_loc[pp]
  }

  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.04)

}
"

FP_hetero_binary <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha) 
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)


  # site-specific intercepts
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- beta0_loc[pp]
  }
  
  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

  ### g-computation in target site
  for (i in 1:N_target) {
    eta[i] <- theta[shared_p + target_site] + inprod(theta[2:beta_p], X_target[i, 2:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


### complete data models ==============
complete_continuous <- "
model{
  ### Likelihood
  prec_y <- 1/exp(log_sigma2)
  for (i in 1:N) {
    mu[i] <- beta0_hetero[site[i]] + alpha[Nperiod + 1 - time[i]] + inprod(theta[1:beta_p], X[i, 1:beta_p])
    y[i] ~ dnorm(mu[i], prec_y)
  }

  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha)
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  # hetero intercept
  for (pp in 1:n_site) {
    beta0_hetero[pp] <- beta0 + delta[pp]
    delta[pp] ~ dnorm(0, tau_delta[pp])
    tau_delta[pp] <- ifelse(pp == target_site, 1.0E10, 1 / (tau2 * lambda2[pp]))
    lambda2[pp] <- pow(lambda[pp], 2)
    lambda[pp] ~ dt(0, 1, 1) T(0,)
  }
  tau2 <- pow(tau, 2)
  tau ~ dt(0, 1, 1) T(0,)

  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.04)

}
"


complete_binary <- "
model{
  ### Likelihood
  for (i in 1:N) {
    eta[i] <- beta0_hetero[site[i]] + alpha[Nperiod + 1 - time[i]] + inprod(theta[1:beta_p], X[i, 1:beta_p])
    p[i] <- ilogit(eta[i])
    y[i] ~ dbern(p[i])
  }

  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha)
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  # hetero intercept
  for (pp in 1:n_site) {
    beta0_hetero[pp] <- beta0 + delta[pp]
    delta[pp] ~ dnorm(0, tau_delta[pp])
    tau_delta[pp] <- ifelse(pp == target_site, 1.0E10, 1 / (tau2 * lambda2[pp]))
    lambda2[pp] <- pow(lambda[pp], 2)
    lambda[pp] ~ dt(0, 1, 1) T(0,)
  }
  tau2 <- pow(tau, 2)
  tau ~ dt(0, 1, 1) T(0,)

  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  
  ### g-computation in target site
  for (i in 1:N_target) {
    eta_g[i] <- beta0_hetero[target_site] + inprod(theta[2:beta_p], X_target[i, 2:beta_p])
    p_ctrl[i] <- ilogit(eta_g[i])
    p_trt[i] <- ilogit(eta_g[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


### time machine models ==============
TM_continuous <- "
model{
  ### Likelihood
  prec_y <- 1/exp(log_sigma2)
  for (i in 1:N) {
    mu[i] <- beta0 + alpha[Nperiod + 1 - time[i]] + inprod(theta[1:beta_p], X[i, 1:beta_p])
    y[i] ~ dnorm(mu[i], prec_y)
  }

  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha)
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.04)

}
"


TM_binary <- "
model{
  ### Likelihood
  for (i in 1:N) {
    eta[i] <- beta0 + alpha[Nperiod + 1 - time[i]] + inprod(theta[1:beta_p], X[i, 1:beta_p])
    p[i] <- ilogit(eta[i])
    y[i] ~ dbern(p[i])
  }

  alpha[1] = 0
  alpha[2] ~ dnorm(0, tau_alpha)
  for(kk in 3:Nperiod) {
      alpha[kk] ~ dnorm(2*alpha[kk - 1] - alpha[kk - 2], tau_alpha)
  }
  tau_alpha <- 1/(sd_alpha)^2
  sd_alpha ~ dt(0, 5^(-2), 1)T(0,)

  ### Prior
  beta0 ~ dnorm(0, 0.01)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  
  ### g-computation in target site
  for (i in 1:N_target) {
    eta_g[i] <- beta0 + inprod(theta[2:beta_p], X_target[i, 2:beta_p])
    p_ctrl[i] <- ilogit(eta_g[i])
    p_trt[i] <- ilogit(eta_g[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"
