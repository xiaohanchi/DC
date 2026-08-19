
# Federated models ==============
## Continous ==============
### All Laplace version ==============
FP_continuous <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])
  
  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }

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
    theta[(shared_p + n_sigma + pp)] <- beta0 + delta[pp]
    delta[pp] ~ dnorm(0, tau_delta[pp])
    tau_delta[pp] <- ifelse(pp == target_site, 1.0E10, 1 / (tau2 * lambda2[pp]))
    lambda2[pp] <- pow(lambda[pp], 2)
    lambda[pp] ~ dt(0, 1, 1) T(0,)
  }
  tau2 <- pow(tau, 2)
  tau ~ dt(0, 1, 1) T(0,)
  
  ### Prior (unified; formerly intercept/noise used 0.01/100)
  beta0 ~ dnorm(0, 0.01)
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

}
"

# borrow on intercept (horseshoe); independent time effects (no RW)
FP_continuous_indepTime <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])
  
  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

  # hetero intercept
  for (pp in 1:n_site) {
    theta[(shared_p + n_sigma + pp)] <- beta0 + delta[pp]
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

}
"

FP_hetero_continuous <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }
  
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

  # site-specific intercepts (unified; formerly 0.01/100)
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    # beta0_loc[pp] ~ dnorm(0, 0.01/100)
    theta[(shared_p + n_sigma + pp)] <- beta0_loc[pp]
  }

  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

}
"

# no borrow on intercept; independent time effects (no RW)
FP_hetero_continuous_indepTime <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }
  
  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

  # site-specific intercepts
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + n_sigma + pp)] <- beta0_loc[pp]
  }

  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

}
"

### Target Exact version ==============
FP_continuous_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact (individual normal)
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  log_s2_exact <- log_sigma2_vec[equals(n_sigma, 1) * 1 + (1 - equals(n_sigma, 1)) * target_site]
  for (i in 1:N_exact) {
    mu_exact[i] <- theta[shared_p + n_sigma + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dnorm(mu_exact[i], exp(-log_s2_exact))
  }

  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }

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
    theta[(shared_p + n_sigma + pp)] <- beta0 + delta[pp]
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

}
"
FP_continuous_indepTime_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  log_s2_exact <- log_sigma2_vec[equals(n_sigma, 1) * 1 + (1 - equals(n_sigma, 1)) * target_site]
  for (i in 1:N_exact) {
    mu_exact[i] <- theta[shared_p + n_sigma + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dnorm(mu_exact[i], exp(-log_s2_exact))
  }
  
  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

  # hetero intercept
  for (pp in 1:n_site) {
    theta[(shared_p + n_sigma + pp)] <- beta0 + delta[pp]
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

}
"

FP_hetero_continuous_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  log_s2_exact <- log_sigma2_vec[equals(n_sigma, 1) * 1 + (1 - equals(n_sigma, 1)) * target_site]
  for (i in 1:N_exact) {
    mu_exact[i] <- theta[shared_p + n_sigma + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dnorm(mu_exact[i], exp(-log_s2_exact))
  }

  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }
  
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
    theta[(shared_p + n_sigma + pp)] <- beta0_loc[pp]
  }

  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

}
"



FP_hetero_continuous_indepTime_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  log_s2_exact <- log_sigma2_vec[equals(n_sigma, 1) * 1 + (1 - equals(n_sigma, 1)) * target_site]
  for (i in 1:N_exact) {
    mu_exact[i] <- theta[shared_p + n_sigma + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dnorm(mu_exact[i], exp(-log_s2_exact))
  }

  for (pp in 1:n_sigma) {
    log_sigma2_vec[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- log_sigma2_vec[pp]
  }
  
  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

  # site-specific intercepts
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + n_sigma + pp)] <- beta0_loc[pp]
  }

  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

}
"

## Binary ==============
### All Laplace version ==============

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
  
  ### Prior (unified; formerly intercept used 0.01/100)
  beta0 ~ dnorm(0, 0.01)
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"

FP_binary_indepTime <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

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

  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

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
    # beta0_loc[pp] ~ dnorm(0, 0.01/100)
    theta[(shared_p + pp)] <- beta0_loc[pp]
  }
  
  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"

FP_hetero_binary_indepTime <- "
model{
  ### Likelihood
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

  # site-specific intercepts
  for (pp in 1:n_site) {
    beta0_loc[pp] ~ dnorm(0, 0.01)
    theta[(shared_p + pp)] <- beta0_loc[pp]
  }
  
  ### Prior
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }

  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


### Target Exact version ==============

FP_binary_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact (individual Bernoulli)
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (i in 1:N_exact) {
    logit_exact[i] <- theta[shared_p + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dbern(ilogit(logit_exact[i]))
  }

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
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"

FP_binary_indepTime_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (i in 1:N_exact) {
    logit_exact[i] <- theta[shared_p + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dbern(ilogit(logit_exact[i]))
  }

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

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
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"

FP_hetero_binary_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (i in 1:N_exact) {
    logit_exact[i] <- theta[shared_p + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dbern(ilogit(logit_exact[i]))
  }

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
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"

FP_hetero_binary_indepTime_exact <- "
model{
  ### Likelihood: non-target Laplace + target exact
  y_laplace[1:n_p] ~ dmnorm(theta[], invSigma[, ])

  for (i in 1:N_exact) {
    logit_exact[i] <- theta[shared_p + target_site] + inprod(X_exact[i, 1:shared_p], theta[1:shared_p])
    y_exact[i] ~ dbern(ilogit(logit_exact[i]))
  }

  for (pp in 1:(Nperiod-1)) {
    theta[(beta_p + pp)] <- alpha[Nperiod - (pp - 1)]
  }
  alpha[1] = 0
  for(kk in 2:Nperiod) {
      alpha[kk] ~ dnorm(0, 0.01)
  }

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
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta[i] <- theta[shared_p + target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta[i])
    p_trt[i] <- ilogit(eta[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


# complete data models ==============
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

  ### Prior (unified; formerly intercept/noise used 0.01/100)
  beta0 ~ dnorm(0, 0.01)
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.01)
  # log_sigma2 ~ dnorm(0, 0.01/100)

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
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  
  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta_g[i] <- beta0_hetero[target_site] + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta_g[i])
    p_trt[i] <- ilogit(eta_g[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


# unadjusted: aggregate by (site, time, trt)
complete_binary_agg <- "
model{
  ### Likelihood
  for (s in 1:N_strata) {
    eta[s] <- beta0_hetero[site[s]] + alpha[Nperiod + 1 - time[s]] + theta[1] * trt[s]
    p[s] <- ilogit(eta[s])
    y[s] ~ dbin(p[s], n[s])
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
  # beta0 ~ dnorm(0, 0.01/100)
  theta[1] ~ dnorm(0, 0.01)

  phat_ctrl <- ilogit(beta0_hetero[target_site])
  phat_trt <- ilogit(beta0_hetero[target_site] + theta[1])

}
"


# time machine models ==============
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
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  log_sigma2 ~ dnorm(0, 0.01)
  # log_sigma2 ~ dnorm(0, 0.01/100)

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
  # beta0 ~ dnorm(0, 0.01/100)
  for (pp in 1:beta_p){
    theta[pp] ~ dnorm(0, 0.01)
  }
  
  ### g-computation in target site (loop avoids invalid 2:beta_p when beta_p = 1)
  for (i in 1:N_target) {
    for (j in 1:beta_p) {
      xcontrib[i, j] <- step(j - 1.5) * theta[j] * X_target[i, j]
    }
    eta_g[i] <- beta0 + sum(xcontrib[i, 1:beta_p])
    p_ctrl[i] <- ilogit(eta_g[i])
    p_trt[i] <- ilogit(eta_g[i] + theta[1])
  }
  phat_ctrl <- sum(p_ctrl[1:N_target]) / N_target
  phat_trt <- sum(p_trt[1:N_target]) / N_target

}
"


# unadjusted: aggregate by (time, trt)
TM_binary_agg <- "
model{
  ### Likelihood
  for (s in 1:N_strata) {
    eta[s] <- beta0 + alpha[Nperiod + 1 - time[s]] + theta[1] * trt[s]
    p[s] <- ilogit(eta[s])
    y[s] ~ dbin(p[s], n[s])
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
  # beta0 ~ dnorm(0, 0.01/100)
  theta[1] ~ dnorm(0, 0.01)

  phat_ctrl <- ilogit(beta0)
  phat_trt <- ilogit(beta0 + theta[1])

}
"
