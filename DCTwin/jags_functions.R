#### MAP.model 

# half-t with scale=10 and df = 2
# dt(0, 0.01, 2)T(0,)
# half-Cauchy is a special case of the t distribution, with df = 1
# scale = 10
# dt(0, 0.01, 1)T(0,)

unadjMAP.normal <- "
model{
  ### RCT
  for (ii in 1:N_RCT) {
    y_rct[ii] ~ dnorm(mu_rct[ii], tau_rct)
    mu_rct[ii] <- mu_ctrl + beta_trt*treatment[ii]
  }
  
  mu_ctrl <- wt * dist1 + (1 - wt) * dist0
  wt ~ dbern(w0)
  dist0 ~ dnorm(0, tau_dist0)
  dist1 ~ dnorm(ybar_syn, tau_dist1)
  tau_dist0 <- 1/var_dist0
  tau_dist1 <- 1/var0
  
  ### Prior
  beta_trt ~ dnorm(0, 0.01)
  tau_rct <- 1/(sd_rct)^2
  sd_rct ~ dt(0, 5^(-2), 1)T(0,)
  
  ### Prediction
  for (ii in 1:N_RCT) {
    y_rct_pred[1, ii] <- mu_ctrl 
    y_rct_pred[2, ii] <- mu_ctrl + beta_trt
  }
}
"

unadjMAP.binary <- "
model{
  ### RCT
  for (ii in 1:N_RCT) {
    y_rct[ii] ~ dbern(p_rct[ii])
    p_rct[ii] <- ilogit(mu_rct[ii])
    mu_rct[ii] <- mu_ctrl + beta_trt*treatment[ii]
  }
  
  mu_ctrl <- wt * dist1 + (1 - wt) * dist0
  wt ~ dbern(w0)
  dist0 ~ dnorm(0, tau_dist0)
  dist1 ~ dnorm(ybar_syn, tau_dist1)
  tau_dist0 <- 1/var_dist0
  tau_dist1 <- 1/var0
  
  ### Prior
  beta_trt ~ dnorm(0, 0.01)
  
  ### Prediction
  for (ii in 1:N_RCT) {
    p_rct_pred[1, ii] <- ilogit(mu_ctrl) 
    p_rct_pred[2, ii] <- ilogit(mu_ctrl + beta_trt)
  }
}
"



adjMAP.normal <- "
model{
  for (pp in 1:P){
    X_mean[pp] <- mean(X[, pp])
  }
  ### RCT
  for (ii in 1:N_RCT) {
    y_rct[ii] ~ dnorm(mu_rct[ii], tau_rct)
    mu_rct[ii] <- mu_ctrl + beta_trt*treatment[ii] + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P]))
  }
  
  mu_ctrl <- wt * dist1 + (1 - wt) * dist0
  wt ~ dbern(w0)
  dist0 ~ dnorm(0, tau_dist0)
  dist1 ~ dnorm(ybar_syn, tau_dist1)
  tau_dist0 <- 1/var_dist0
  tau_dist1 <- 1/var0
  
  ### Prior
  beta_trt ~ dnorm(0, 0.01)
  for (pp in 1:P){
    beta[pp] ~ dnorm(0, 0.01)
  }
  tau_rct <- 1/(sd_rct)^2
  sd_rct ~ dt(0, 5^(-2), 1)T(0,)
  
  ### Prediction
  for (ii in 1:N_RCT) {
    y_rct_pred[1, ii] <- mu_ctrl + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P]))
    y_rct_pred[2, ii] <- mu_ctrl + beta_trt + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P]))
  }
}
"

DCTwin.binary <- "
model{
  ### Syn Data
  for (jj in 1:Ngroup) {
    ysum_syn[jj] ~ dbin(p_syn[jj], n_dc)
    p_syn[jj] <- ilogit(mu_syn[jj])
    mu_syn[jj] ~ dnorm(mu, tau0_syn)
  }
  
  ### RCT: sufficient likelihood
  muhat_ctrl ~ dnorm(mu_ctrl, tau_ctrl)
  muhat_trt ~ dnorm(mu_trt, tau_trt)
  
  mu_ctrl <- wt * dist1 + (1 - wt) * dist0
  wt ~ dbern(w0)
  dist0 ~ dnorm(0, tau_dist0)
  dist1 ~ dnorm(mu, tau_dist1)
  tau_dist0 <- 1/var_dist0
  tau_dist1 <- 1/var_dist1
  var_dist1 <- var0 + (sd0_syn)^2
  
  ### Prior
  mu ~ dnorm(0, 0.01)
  tau0_syn <- 1/(sd0_syn)^2
  sd0_syn ~ prior_to_be_defined
  mu_trt ~ dnorm(0, 0.01)
  
  ### Prediction
  p_est[1] <- ilogit(mu_ctrl)
  p_est[2] <- ilogit(mu_trt)
  
}
"

logistic_RCT <- "
model{
  for (pp in 1:P){
      X_mean[pp] <- mean(X[, pp])
  }
    
  for (ii in 1:N_RCT) {
      y_rct[ii] ~ dbern(p_rct[ii])
      p_rct[ii] <- ilogit(mu_rct[ii])
      mu_rct[ii] <- mu_ctrl + beta_trt*treatment[ii] + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P]))
  }
  
  mu_ctrl ~ dnorm(0, 0.01)
  beta_trt ~ dnorm(0, 0.01)
  for (pp in 1:P){
    beta[pp] ~ dnorm(0, 0.01)
  }
    
  ### Prediction
  for (ii in 1:N_RCT) {
    p_ctrl_pred[ii] <- ilogit(mu_ctrl + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P])))
    p_trt_pred[ii] <- ilogit(mu_ctrl + beta_trt + sum(beta[1:P] * (X[ii, 1:P] - X_mean[1:P])))
  }
}
"
