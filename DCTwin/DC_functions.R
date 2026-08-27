#### data preparation ===============

logit <- function(p) log(p / (1 - p))
expit <- function(x) exp(x) / (1 + exp(x))

wt.func <- function(x, ref.stat = qnorm(1-0.025), rho = 10, b = 2, type = 2){
  if(type == 1) {
    a <- rho/ref.stat^b
    output <- 1/exp(a * abs(x)^b)
  } else if (type == 2) {
    a <- b/rho
    output <- 1/(1 + exp(a*abs(x) - b))
  } else if (type == 3) {
    a <- rho/(1 - rho)
    num <- a * abs(x)^b
    den <- num + (1 - abs(x))^b
    output <- num / den
  }
  return(output)
}

#### simulation func =============
get.Gphi.weight <- function(data, Gphi.model, n_simu = 1000) {
  X_v <- data.matrix(data)
  disc_col <- which(apply(X_v, 2, function(z) length(unique(z)) <= 10L))
  if (length(disc_col)) {
    X_v[, disc_col] <- sweep(
      X_v[, disc_col, drop = FALSE], 2, apply(X_v[, disc_col, drop = FALSE], 2, min)
    )
  }
  X_tilde <- reticulate::py_to_r(
    Gphi.model$generate_synthetic_data(
      n_samples = as.integer(n_simu), t = 1, n_permutations = 2L
    )$detach()$cpu()$numpy()
  )
  n_v <- nrow(X_v)
  m_G <- nrow(X_tilde)
  kappa_G <- 0.5
  cross_d <- sqrt(pmax(
    outer(rowSums(X_v^2), rowSums(X_tilde^2), "+") - 2 * tcrossprod(X_v, X_tilde),
    0
  ))
  E_hat <- (2 / (n_v * m_G)) * sum(cross_d) -
    (2 / n_v^2) * sum(dist(X_v)) -
    (2 / m_G^2) * sum(dist(X_tilde))
  D_hat <- E_hat / mean(dist(rbind(X_v, X_tilde)))
  n_eff <- n_v * m_G / (n_v + m_G)
  exp(-n_eff^kappa_G * D_hat)
}

digital.control <- function(rwd.data, exp.all, EHR.data, RCT.data, 
                            synctrl.n, seed, Gphi.type = c(1, 2, 3)) {
  # synctrl.n: sample size of synthetic controls in stage 1 (e.g., 100)
  # Gphi.type = 1 for using RCT + RWD data in structure learning; Gphi.type = 2 for using RCT data only
  set.seed(seed)
  
  rwd.data <- rwd.data %>% select(-c(treatment))
  exp.data <- RCT.data %>% filter(treatment == 1)
  ctrl.data <- RCT.data %>% filter(treatment == 0)
  EHR.data <- EHR.data %>% select(-c(treatment, S, Y))
  
  ### TabPFN prediction for DC: using only RWD
  pred.model <- tab_pfn(
    Y ~ ., data = rwd.data, version = "v3", training_set_limit = Inf,
    control = control_tab_pfn(device = "auto", n_preprocessing_jobs = 1L, ignore_pretraining_limits = TRUE)
  )
  y.pred.dt <- predict(
    pred.model, new_data = (ctrl.data %>% select(-c(treatment, S, Y)))
  )$.pred
  
  ### TabPFN for X covariates (same training X as BN parameter learning)
  if (Gphi.type == 1) {
    Gphi.data <- RCT.data %>% select(-c(treatment, S, Y))
  } else if (Gphi.type == 2) {
    Gphi.data <- EHR.data
  } else if (Gphi.type == 3) {
    Gphi.data <- RCT.data %>% select(-c(treatment, S, Y))
  }
  coltype <- sapply(1:ncol(Gphi.data), function(r) class(data.frame(Gphi.data)[, r]))
  rwd.data <- rwd.data %>% data.frame()
  tabpfn_py <- reticulate::import("tabpfn")
  Gphi.model <- reticulate::import("tabpfn_extensions.unsupervised")$TabPFNUnsupervisedModel(
    tabpfn_clf = tabpfn_py$TabPFNClassifier(device = "auto"),
    tabpfn_reg = tabpfn_py$TabPFNRegressor(device = "auto")
  )
  X <- data.matrix(Gphi.data)
  disc_col <- which(apply(X, 2, function(z) length(unique(z)) <= 10L))
  if (length(disc_col)) {
    X[, disc_col] <- sweep(X[, disc_col, drop = FALSE], 2, apply(X[, disc_col, drop = FALSE], 2, min))
    Gphi.model$set_categorical_features(as.list(as.integer(disc_col - 1L)))
  }
  Gphi.model$fit(X)
    
  ### generate X
  Gphi.ctrl2 <- as_tibble(
    reticulate::py_to_r(
      Gphi.model$generate_synthetic_data(n_samples = as.integer(synctrl.n), t = 1, n_permutations = 2L)$detach()$cpu()$numpy()
    ),
    .name_repair = ~ names(Gphi.data)
  )
  if (sum(coltype != "numeric") == 3) {
    Gphi.ctrl2  <- Gphi.ctrl2 %>%
      mutate(X5 = factor(X5), X6 = factor(X6), X7 = factor(X7))
  } else if (sum(coltype != "numeric") == 2) {
    Gphi.ctrl2 <- Gphi.ctrl2 %>% mutate(X6 = factor(X6), X7 = factor(X7))
  }

  full_out <- pred.model$fit$predict(
    hardhat::forge(as.data.frame(Gphi.ctrl2), pred.model$blueprint)$predictors,
    output_type = "full"
  )
  y.pred.syn1 <- matrix(as.numeric(full_out$criterion$sample(full_out$logits)$cpu()$numpy()), ncol = 1)
  
  return(list(
    Gphi.model = Gphi.model, 
    Gphi.ctrl2 = Gphi.ctrl2,
    y.pred.dt = y.pred.dt, 
    y.pred.syn1 = y.pred.syn1
  ))
}

get.DC.ATE <- function(jagsmodel, RCT.data, type, cross.fit = FALSE){
  var.name <- ifelse(type == 1, "y_rct_pred", "p_rct_pred")
  y.pred0.sample <- as.mcmc.list(
    jagsmodel,
    sapply(1:nrow(RCT.data), function(r) paste0(var.name, "[1,", r, "]"))
  ) %>% as.matrix()
  y.pred1.sample <- as.mcmc.list(
    jagsmodel,
    sapply(1:nrow(RCT.data), function(r) paste0(var.name, "[2,", r, "]"))
  ) %>% as.matrix()
  # dirichlet weight: Bayesian bootstrap
  weight <- rdirichlet(1, rep(1, nrow(RCT.data)))
  y.pred0.dist <- as.vector(y.pred0.sample %*% t(weight))
  y.pred1.dist <- as.vector(y.pred1.sample %*% t(weight))
  
  ATE.mcmc <- y.pred1.dist - y.pred0.dist
  if (!cross.fit) {
    results <- c(
      prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)), 
      ATE = mean(ATE.mcmc)
    )
  } else if (cross.fit) {
    results <- list(
      ATE.mcmc = ATE.mcmc
    )
  }
  return(results)
}

#### MAP function ======
MAP.func <- function(rawRWD, RCT.data, true.ctrl.s1, exp.all, 
                     var0.ess, wt.rho.x, wt.b.x, wt.rho.y, wt.b.y, wt.type, w0.val, 
                     RCT.data.trans, Gphi.model, DC.groups, y.pred.dt, y.pred.syn1, 
                     methods = c("full", "selected"), seed) {
  # var0.ess: sigma0^2 ESS value for the non-informative part in the mixture prior
  # w0.val: = 2 or 3 for varying weighting function
  set.seed(seed)
  methods <- match.arg(methods)
  runjags.options(silent.jags = TRUE, silent.runjags = TRUE, inits.warning = FALSE)
  
  RWD <- rawRWD %>% select(-c(treatment))
  true.ctrl.s1 <- true.ctrl.s1 %>% select(-c(label))
  ctrl.data <- RCT.data %>% filter(treatment == 0) %>% select(., -c(treatment, S))
  exp.all <- exp.all %>% select(-c(label))
  ctrl.all <- bind_rows(true.ctrl.s1, ctrl.data)
  
  outcome.type <- ifelse(length(unique(RWD$Y)) == 2, 2, 1) # 2 for binary; 1 for continuous
  
  ######################################## DC ########################################
  
  # get weights from ML models: self-validation
  wt.para <- c()
  t_stat <- t.test(y.pred.dt, ctrl.data$Y, paired = TRUE)$statistic
  pval <- 2*(1 - pnorm(abs(t_stat) * length(y.pred.dt)^(-0.01)))
  wt.para[1] <- wt.func(x = pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
  
  Gphi.pval <- get.Gphi.weight(
    data = RCT.data %>% select(-c(treatment, S, Y)), Gphi.model = Gphi.model
  )
  wt.para[2] <- wt.func(x = Gphi.pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
  
  wt.para <- pmax(wt.para, 1e-12)
  
  if(outcome.type == 1) {
    ### DC_unadj_v2: std
    jagsdata <- list(
      ybar_syn = as.numeric(colMeans(y.pred.syn1)),
      N_RCT = nrow(RCT.data),
      var0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/nrow(ctrl.data), 
      var_dist0 =  max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/var0.ess,
      y_rct = RCT.data$Y, 
      treatment = RCT.data$treatment, 
      w0 = case_when(
        w0.val == -1 ~ min(wt.para),
        w0.val == 2 ~ wt.para[1],
        w0.val == 3 ~ wt.para[2],
        TRUE ~ w0.val
      )
    )
    jagsmodel <- run.jags(
      model = unadjMAP.normal, 
      monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"), 
      data = jagsdata, n.chains = 4, 
      adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      method = "rjags", plots = FALSE, silent.jags = T,
      inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    )
    
    res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data, type = outcome.type)
    prob.DCunadj2 <- res.tmp["prob"]
    ATE.DCunadj2 <- res.tmp["ATE"]
    
    ### ATE: DC_adj_v1
    jagsdata <- list(
      ybar_syn = as.numeric(colMeans(y.pred.syn1)),
      N_RCT = nrow(RCT.data),
      P = sum(grepl("X", colnames(RCT.data))), 
      var0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/nrow(ctrl.data), 
      var_dist0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/var0.ess,
      y_rct = RCT.data$Y, 
      treatment = RCT.data$treatment, 
      w0 = case_when(
        w0.val == -1 ~ min(wt.para),
        w0.val == 2 ~ wt.para[1],
        w0.val == 3 ~ wt.para[2],
        TRUE ~ w0.val
      ),
      X = (apply(as.matrix(RCT.data), c(1, 2), as.numeric)[, grep("X", names(RCT.data)), drop = FALSE])
    )
    
    jagsmodel <- run.jags(
      model = adjMAP.normal, 
      monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"), 
      data = jagsdata, n.chains = 4, 
      adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      method = "rjags", plots = FALSE, silent.jags = T, 
      inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    )
    
    res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data, type = outcome.type)
    prob.DCadj1 <- res.tmp["prob"]
    ATE.DCadj1 <- res.tmp["ATE"]
    
    ### ATE: DC_twin: only use m(x) in covariate adjustment
    prog.model <- tab_pfn(
      Y ~ ., data = RWD, version = "v3", training_set_limit = Inf,
      control = control_tab_pfn(device = "auto", n_preprocessing_jobs = 1L, ignore_pretraining_limits = TRUE)
    )
    prog.score <- predict(
      prog.model, new_data = (RCT.data %>% select(-c(S, treatment, Y)))
    )$.pred
    
    jagsdata <- list(
      ybar_syn = as.numeric(colMeans(y.pred.syn1)),
      N_RCT = nrow(RCT.data),
      P = 1, 
      var0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/nrow(ctrl.data), 
      var_dist0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/var0.ess,
      y_rct = RCT.data$Y, 
      treatment = RCT.data$treatment, 
      w0 = case_when(
        w0.val == -1 ~ min(wt.para),
        w0.val == 2 ~ wt.para[1],
        w0.val == 3 ~ wt.para[2],
        TRUE ~ w0.val
      ),
      X = cbind(prog.score)
    )
    
    jagsmodel <- run.jags(
      model = adjMAP.normal, 
      monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"), 
      data = jagsdata, n.chains = 4, 
      adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      method = "rjags", plots = FALSE, silent.jags = T, 
      inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    )
    
    res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data, type = outcome.type)
    prob.DCtwin <- res.tmp["prob"]
    ATE.DCtwin <- res.tmp["ATE"]
    
  } else if (outcome.type == 2) {
    
    # to be updated
    
    jagsdata <- list(
      ybar_syn = logit(as.numeric(pmin(pmax(colMeans(y.pred.syn1), 1e-6), 1 - 1e-6))),
      N_RCT = nrow(RCT.data),
      var0 = 1/(max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y)) * nrow(y.pred.syn1)), 
      var_dist0 = 1/(max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y)) * var0.ess),
      y_rct = RCT.data$Y, 
      treatment = RCT.data$treatment, 
      w0 = case_when(
        w0.val == -1 ~ min(wt.para),
        w0.val == 2 ~ wt.para[1],
        w0.val == 3 ~ wt.para[2],
        TRUE ~ w0.val
      )
    )
    jagsmodel <- run.jags(
      model = unadjMAP.binary, 
      monitor = c("mu_ctrl", "beta_trt", "p_rct_pred"), 
      data = jagsdata, n.chains = 4, 
      adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      method = "rjags", plots = FALSE, silent.jags = T,
      inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    )
    
    res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data, type = outcome.type)
    prob.DCunadj2 <- res.tmp["prob"]
    ATE.DCunadj2 <- res.tmp["ATE"]
    
    ### ATE: DC_adj_v1
    prob.DCadj1 <- NA
    ATE.DCadj1 <- NA
    
    prob.DCtwin <- NA
    ATE.DCtwin <- NA
    
  }
  
  
  ################################### Competing Methods ##################################
  ### ATE: RCT
  ATE.anova <- mean(exp.all$Y) - mean(ctrl.data$Y)
  pval.anova <- t.test(exp.all$Y, ctrl.data$Y, alternative = "two.sided")$p.value
  
  ### ATE: 1:1 RCT
  ATE.anova2 <- mean(exp.all$Y) - mean(ctrl.all$Y)
  pval.anova2 <- t.test(exp.all$Y, ctrl.all$Y, alternative = "two.sided")$p.value
  
  ### ANCOVA
  if(outcome.type == 1) {
    rct.fit <- lm(Y ~ ., data = (RCT.data %>% select(-c(S))))
    ATE.ancova <- rct.fit$coefficients["treatment"]
    pval.ancova <- 2 * (1 - pt(abs(coef(summary(rct.fit))["treatment", 3]), rct.fit$df))
  } else if (outcome.type == 2) {
    rct.fit <- glm(Y ~ ., data = (RCT.data %>% select(-c(S))), family = binomial)
    ATE.ancova <- mean(predict(rct.fit, newdata = (RCT.data %>% mutate(treatment = 1)), type = "response")) - mean(predict(rct.fit, newdata = (RCT.data %>% mutate(treatment = 0)), type = "response"))
    pval.ancova <- coef(summary(rct.fit))["treatment", "Pr(>|z|)"]
  }
  
  if (methods == "full") {
    ### ATE: PSPP with all current data
    pspp.data <- bind_rows(
      (ctrl.data %>% mutate(label = 1, arm = 0)),
      (exp.all %>% mutate(label = 1, arm = 1)),
      (RWD %>% mutate(label = 0, arm = 0))
    )
    
    for (ii in 5:1) {
      pspp.try <- try(ATE.pspp(
        pspp.data = pspp.data, strata.n = ii, borrow.n = nrow(ctrl.data), 
        type = 1, outcome.type = ifelse(outcome.type == 1, "continuous", "binary")
      ), silent = TRUE)
      if (!("try-error" %in% class(pspp.try))) break
    }
    ATE.pspp1 <- pspp.try$ATE
    prob.pspp1 <- pspp.try$Prob
    
    ### ATE: PSMAP with all current data
    for (ii in 5:1) {
      psmap.try <- try(ATE.psmap(
        hist.data = RWD, current.ctrl = ctrl.data, current.trt = exp.all,
        MAP_ESS = nrow(ctrl.data), S = ii, showESS = FALSE, type = 1, outcome.type = outcome.type
      ), silent = TRUE)
      if (!("try-error" %in% class(psmap.try))) break
    }
    ATE.psmap1 <- psmap.try["ATE"]
    prob.psmap1 <- psmap.try["prob"]
    
    ### ATE: PW-MEM: borrow n_ctrl patients (because of PW of control patients)
    tmp.pwmem <- ATE.pwmem(
      RWD = RWD, ctrl.data = ctrl.data, exp.all = exp.all, outcome.type = outcome.type
      )
    ATE.pwmem <- tmp.pwmem$ATE
    prob.pwmem <- tmp.pwmem$Prob
    
    ### ATE: PS-SAM (CSD = 0.4) 
    # 3: 2-fold cross fitting
    idx2_0 <- sample(1:nrow(RWD), (nrow(RWD)/2), replace = FALSE)
    idx3_0 <- sample(1:nrow(ctrl.data), (nrow(ctrl.data)/2), replace = FALSE)
    idx4_0 <- sample(1:nrow(exp.all), (nrow(exp.all)/2), replace = FALSE)
    for (fold in 1:2) {
      if(fold == 1) {
        idx2 <- idx2_0
        idx3 <- idx3_0
        idx4 <- idx4_0
      } else {
        idx2 <- setdiff(1:nrow(RWD), idx2_0)
        idx3 <- setdiff(1:nrow(ctrl.data), idx3_0)
        idx4 <- setdiff(1:nrow(exp.all), idx4_0)
      }
      RCT.subset <- bind_rows(
        (ctrl.data[-idx3, ] %>% mutate(label = 1, arm = 0)),
        (exp.all[-idx4, ] %>% mutate(label = 1, arm = 1))
      ) %>% dplyr::sample_n(nrow(ctrl.data)/2, replace = TRUE)
      ps.data <- bind_rows(RCT.subset, (RWD %>% mutate(label = 0, arm = 0)))
      tmp.pssam <- ATE.pssam.cf(
        ps.data = ps.data, ctrl.data.valid = ctrl.data[-idx3, ], ctrl.data.est = ctrl.data[idx3, ], 
        exp.data = exp.all[idx4, ], sigma = sd(RWD$Y), eff.size = ifelse(outcome.type == 1, 0.4, 0.15), 
        outcome.type = outcome.type
      )
      if (fold == 1) {
        ATE.mcmc <- tmp.pssam$ATE.mcmc
      } else {
        ATE.mcmc <- cbind(ATE.mcmc, tmp.pssam$ATE.mcmc)
      }
    }
    ATE.mcmc.pool <- rowMeans(ATE.mcmc)
    prob.pssam3 = 2 * min(mean(ATE.mcmc.pool < 0), mean(ATE.mcmc.pool > 0))
    ATE.pssam3 = mean(ATE.mcmc.pool)
    
    
    # 4: borrow n_dc patients
    RCT.subset <- bind_rows(
      (ctrl.data %>% mutate(label = 1, arm = 0)),
      (exp.all %>% mutate(label = 1, arm = 1))
    ) %>% dplyr::sample_n(nrow(ctrl.data), replace = TRUE)
    
    ps.data <- bind_rows(RCT.subset, (RWD %>% mutate(label = 0, arm = 0)))
    tmp.pssam <- ATE.pssam(
      ps.data = ps.data, ctrl.data = ctrl.data, exp.data = exp.all, 
      sigma = sd(RWD$Y), eff.size = ifelse(outcome.type == 1, 0.4, 0.15), 
      outcome.type = outcome.type
    )
    ATE.pssam4 <- tmp.pssam$ATE
    prob.pssam4 <- tmp.pssam$Prob
    
    # ### original PROCOVA model: trained in RWD
    # prog.model <- bart2(
    #   Y ~ ., data = RWD,  n.trees = 150, n.samples = 2500, n.chains = 4, 
    #   keepTrees = TRUE, combineChains = T, n.threads = 1, verbose = FALSE, seed = 233
    # )
    # prog.score <- predict(
    #   prog.model, select((RCT.data %>% select(-c(S))), -c(treatment, Y))
    # ) %>% colMeans()
    # ### test PROCOVA model: trained in RCT
    # prog.model.rct <- bart2(
    #   Y ~ ., data = ctrl.data, n.trees = 150, n.samples = 2500, n.chains = 4, 
    #   keepTrees = TRUE, combineChains = T, n.threads = 1, verbose = FALSE, seed = 233
    # )
    # ### test PROCOVA model: linear prog model
    # prog.model.lr <- lm(Y ~ ., data = RWD)
    # 
    ### ATE: PROCOVA
    procova.data <- RCT.data %>% select(-c(S)) %>% mutate(prog.score = prog.score)
    for (ii in 1:(ncol(procova.data) - 3)) {
      if (class(procova.data[[ii]]) == "numeric") {
        procova.data[ii] <- procova.data[ii] - mean(pull(procova.data[ii]))
      } else if (class(procova.data[[ii]]) == "factor") {
        numeric_values <- as.numeric(as.character(procova.data[[ii]]))
        procova.data[[ii]] <- numeric_values - mean(numeric_values)
      }
    }
    procova.data$prog.score <- procova.data$prog.score - mean(procova.data$prog.score)
    if(outcome.type == 1) {
      procova.fit <- lm(Y ~ ., data = procova.data)
      ATE.procova <- procova.fit$coefficients["treatment"]
      pval.procova <- coef(summary(procova.fit))["treatment", "Pr(>|t|)"]
    } else if (outcome.type == 2) {
      procova.fit <- glm(Y ~ ., data = procova.data, family = binomial)
      ATE.procova <- mean(predict(procova.fit, newdata = (procova.data %>% mutate(treatment = 1)), type = "response")) - mean(predict(procova.fit, newdata = (procova.data %>% mutate(treatment = 0)), type = "response"))
      pval.procova <- coef(summary(procova.fit))["treatment", "Pr(>|z|)"]
    }
    
    
    # 
    # ### ATE: PROCOVA.rct (only use RCT data)
    # procova.data <- RCT.data %>% select(-c(S))
    # procova.data$prog.score <- predict(prog.model.rct, select(procova.data, -c(treatment, Y))) %>% colMeans()
    # procova.data <- procova.data %>% mutate(treatment = factor(treatment))
    # for (ii in 1:(ncol(procova.data) - 3)) {
    #   if (class(procova.data[[ii]]) == "numeric") {
    #     procova.data[ii] <- procova.data[ii] - mean(pull(procova.data[ii]))
    #   }
    # }
    # procova.data$prog.score <- procova.data$prog.score - mean(procova.data$prog.score)
    # procova.rct.fit <- lm(Y ~ ., data = procova.data)
    # ATE.procova.rct <- procova.rct.fit$coefficients["treatment1"]
    # pval.procova.rct <- pt(coef(summary(procova.rct.fit))["treatment1", 3], procova.rct.fit$df, lower = FALSE)
    # 
    # ### ATE: PROCOVA.lr (lr in RWD)
    # procova.data <- RCT.data %>% select(-c(S))
    # procova.data$prog.score <- predict(prog.model.lr, select(procova.data, -c(treatment, Y)))
    # procova.data <- procova.data %>% mutate(treatment = factor(treatment))
    # for (ii in 1:(ncol(procova.data) - 3)) {
    #   if (class(procova.data[[ii]]) == "numeric") {
    #     procova.data[ii] <- procova.data[ii] - mean(pull(procova.data[ii]))
    #   }
    # }
    # procova.data$prog.score <- procova.data$prog.score - mean(procova.data$prog.score)
    # procova.lr.fit <- lm(Y ~ ., data = procova.data)
    # ATE.procova.lr <- procova.lr.fit$coefficients["treatment1"]
    # pval.procova.lr <- pt(coef(summary(procova.lr.fit))["treatment1", 3], procova.lr.fit$df, lower = FALSE)
    
    ### ATE: Semi-Synthetic control (match 200 -> randomly select 100)
    # for homo population, ATT = ATE
    semi.matchdata <- bind_rows(
      (exp.all %>% mutate(label = 1)), (RWD %>% mutate(label = 0))
    )
    #### ratio = 1:1
    match.it.semi <- matchit(
      label ~ . - Y, data = semi.matchdata, method = "nearest",
      caliper = 0.2, std.caliper = TRUE, ratio = 1, tol = 1e-10
    )
    matchdata.semi <- match.data(match.it.semi)[1:ncol(semi.matchdata)]
    match.ctrl <- filter(matchdata.semi, label == 0) %>% select(-c(label))
    match.ctrl <- match.ctrl[sample(c(1:nrow(match.ctrl)), size = (nrow(exp.all) - nrow(ctrl.data)), replace = FALSE), ]
    ctrl.semiSC <- bind_rows(match.ctrl, ctrl.data)
    ATE.semiSC1 <- mean(exp.all$Y) - mean(ctrl.semiSC$Y)
    pval.semiSC1 <- t.test(exp.all$Y, ctrl.semiSC$Y, alternative = "two.sided")$p.value
    
    # Semi-Synthetic control 2
    semi.matchdata <- bind_rows(
      (ctrl.data %>% mutate(label = 1)), 
      (exp.all %>% mutate(label = 1)), 
      (RWD %>% mutate(label = 0))
    )
    match.it.semi <- matchit(
      label ~ . - Y, data = semi.matchdata, method = "nearest",
      caliper = 0.2, std.caliper = TRUE, ratio = 1, tol = 1e-10
    )
    matchdata.semi <- match.data(match.it.semi)[1:ncol(semi.matchdata)]
    match.ctrl <- filter(matchdata.semi, label == 0) %>% select(-c(label))
    match.ctrl <- match.ctrl[sample(c(1:nrow(match.ctrl)), size = (nrow(exp.all) - nrow(ctrl.data)), replace = FALSE), ]
    ctrl.semiSC <- bind_rows(match.ctrl, ctrl.data)
    ATE.semiSC2 <- mean(exp.all$Y) - mean(ctrl.semiSC$Y)
    pval.semiSC2 <- t.test(exp.all$Y, ctrl.semiSC$Y, alternative = "two.sided")$p.value
    
    
    ### ATE: g-computation (BART) with only RCT trial data
    # bart.rct <- bart2(
    #   Y ~ ., data = (RCT.data %>% select(-c(S))),
    #   n.trees = 150, n.samples = 2500, n.chains = 4, keepTrees = TRUE,
    #   combineChains = T, n.threads = 1, verbose = FALSE, seed = 233
    # )
    # y.pred0.sample <- predict(
    #   bart.rct, (RCT.data %>% select(-c(S)) %>% mutate(treatment = 0) %>% select(-c(Y)))
    # )
    # y.pred1.sample <- predict(
    #   bart.rct, (RCT.data %>% select(-c(S)) %>% mutate(treatment = 1) %>% select(-c(Y)))
    # )
    # 
    # ATE.RCTbart <- mean(colMeans(y.pred1.sample)) - mean(colMeans(y.pred0.sample))
    # pval.RCTbart <- mean(rowMeans(y.pred1.sample) - rowMeans(y.pred0.sample) < 0)
    # 
    # ### ATE: ancova (heteroskedasticity) with RCT & RWD data
    # pooldata <- bind_rows(RCT.data, (rawRWD %>% mutate(S = 1)))
    # lm.fit <- lm(Y ~ ., data = pooldata)
    # robust.se <- vcovHC(lm.fit, type = "HC0")
    # lm.robust <- coeftest(lm.fit, vcov = robust.se)
    # ATE.pool.lmrobust <- coef(lm.robust)["treatment"]
    # pval.pool.lmrobust <- pt(lm.robust["treatment", 3], lm.fit$df, lower = FALSE)
    # 
    # ### ATE: g-computation (BART) with RCT & RWD data
    # bart.model2 <- bart2(
    #   Y ~ ., data = pooldata, n.trees = 150, n.samples = 2500, n.chains = 4, 
    #   keepTrees = TRUE, combineChains = T, n.threads = 1, verbose = FALSE, seed = 233
    # )
    # y.pred0.sample <- predict(
    #   bart.model2, (RCT.data %>% mutate(treatment = 0) %>% select(-c(Y)))
    # )
    # y.pred1.sample <- predict(
    #   bart.model2, (RCT.data %>% mutate(treatment = 1) %>% select(-c(Y)))
    # )
    # ATE.poolbart <- mean(colMeans(y.pred1.sample)) - mean(colMeans(y.pred0.sample))
    # pval.poolbart <- mean(rowMeans(y.pred1.sample) - rowMeans(y.pred0.sample) < 0)
    
    ATE.res <- tibble(
      ATE.semiSC1 = ATE.semiSC1,
      ATE.semiSC2 = ATE.semiSC2,
      # ATE.pw = ATE.pw,
      ATE.pspp1 = ATE.pspp1,
      ATE.psmap1 = ATE.psmap1,
      ATE.pwmem = ATE.pwmem,
      ATE.pssam3 = ATE.pssam3, 
      ATE.pssam4 = ATE.pssam4, 
      # ATE.DCunadj1 = ATE.DCunadj1,
      ATE.DCunadj2 = ATE.DCunadj2,
      ATE.DCadj1 = ATE.DCadj1,
      ATE.DCtwin = ATE.DCtwin, 
      ATE.anova = ATE.anova,
      ATE.anova2 = ATE.anova2,
      # ATE.anova3 = ATE.anova3,
      ATE.ancova = ATE.ancova,
      # ATE.RCTbart = ATE.RCTbart,
      # ATE.pool.lmrobust = ATE.pool.lmrobust,
      # ATE.poolbart = ATE.poolbart,
      ATE.procova = ATE.procova
      # ATE.procova.rct = ATE.procova.rct,
      # ATE.procova.lr = ATE.procova.lr
    )
    
    prob.res <- tibble(
      pval.semiSC1 = pval.semiSC1,
      pval.semiSC2 = pval.semiSC2,
      # pval.pw = pval.pw,
      prob.pspp1 = prob.pspp1,
      prob.psmap1 = prob.psmap1,
      prob.pwmem = prob.pwmem,
      prob.pssam3 = prob.pssam3, 
      prob.pssam4 = prob.pssam4, 
      # prob.DCunadj1 = prob.DCunadj1,
      prob.DCunadj2 = prob.DCunadj2,
      prob.DCadj1 = prob.DCadj1,
      prob.DCtwin = prob.DCtwin,
      pval.anova = pval.anova,
      pval.anova2 = pval.anova2,
      # pval.anova3 = pval.anova3,
      pval.ancova = pval.ancova,
      # pval.RCTbart = pval.RCTbart,
      # pval.pool.lmrobust = pval.pool.lmrobust,
      # pval.poolbart = pval.poolbart,
      pval.procova = pval.procova
      # pval.procova.rct = pval.procova.rct,
      # pval.procova.lr = pval.procova.lr
    )
  } else if (methods == "selected") {
    ATE.res <- tibble(
      ATE.anova = ATE.anova,
      ATE.anova2 = ATE.anova2,
      ATE.ancova = ATE.ancova,
      ATE.DCunadj2 = ATE.DCunadj2,
      ATE.DCadj1 = ATE.DCadj1,
      ATE.DCtwin = ATE.DCtwin,
      pval = pval
    )
    prob.res <- tibble(
      pval.anova = pval.anova,
      pval.anova2 = pval.anova2,
      pval.ancova = pval.ancova,
      prob.DCunadj2 = prob.DCunadj2,
      prob.DCadj1 = prob.DCadj1,
      prob.DCtwin = prob.DCtwin
    )
  }
  
  return(list(ATE = ATE.res, Prob = prob.res))
}



### MAIN function ==================================
MAIN.func <- function(rwd.n, exp.n, EHR.n, synctrl.n, trt.eff, bias.c,
                      scenario, var0.ess, 
                      wt.rho.x, wt.b.x, wt.rho.y, wt.b.y, wt.type, w0.val, 
                      sigma.rwdx = 1, sigma.rwd = 1, 
                      sigma.rctx = 1, sigma.rct = 1, rho.rwd = 0.3, 
                      model.type, bias.type, Gphi.type, outcome.type, seed, rep) {
  # var0.ess: sigma0^2 ESS value for the non-informative part in the mixture prior

  ### Data generating
  tmp.data <- prepare.data(
    rwd.n = rwd.n, exp.n = exp.n, EHR.n = EHR.n, 
    trt.eff = trt.eff, bias.c = bias.c,
    scenario = scenario, 
    sigma.rwdx = sigma.rwdx, sigma.rwd = sigma.rwd, 
    sigma.rctx = sigma.rctx, sigma.rct = sigma.rct, rho.rwd = rho.rwd,
    model.type = model.type, bias.type = bias.type, outcome.type = outcome.type, 
    seed = seed
  )
  rawRWD <- (tmp.data$rawRWD %>% dplyr::select(-c(S)))
  exp.all <- tmp.data$exp.all %>% mutate(label = 1)
  true.ctrl.s1 <- tmp.data$true.ctrl.s1 %>% mutate(label = 1)
  
  res.s1 <- digital.control(
    rwd.data = rawRWD, 
    exp.all = exp.all, 
    EHR.data = tmp.data$EHR.data, 
    RCT.data = tmp.data$RCT.data, 
    synctrl.n = synctrl.n, 
    seed = seed, 
    Gphi.type = Gphi.type
  )
  ### MAP
  res.s2 <- MAP.func(
    rawRWD = rawRWD, 
    RCT.data = tmp.data$RCT.data,
    true.ctrl.s1 = true.ctrl.s1, 
    exp.all = exp.all, 
    var0.ess = var0.ess,
    wt.rho.x = wt.rho.x, 
    wt.b.x = wt.b.x, 
    wt.rho.y = wt.rho.y, 
    wt.b.y = wt.b.y, 
    wt.type = wt.type, 
    w0.val = w0.val, 
    RCT.data.trans = res.s1$RCT.data.trans, 
    Gphi.model = res.s1$Gphi.model, 
    DC.groups = res.s1$Gphi.ctrl2,
    y.pred.dt = res.s1$y.pred.dt, 
    y.pred.syn1 = res.s1$y.pred.syn1, 
    methods = c("full", "selected")[1],
    seed = seed
  )

  output <- list(
    ATE = mutate(res.s2$ATE, Replicate = rep, .before = 1),
    Prob = mutate(res.s2$Prob, Replicate = rep, .before = 1)
  )

  return(output)
}

