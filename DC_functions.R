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

data.transform <- function(x, range = c(), type = "unif_norm") {
  # range: real unif limits when transforming data from std normal to (unscaled) uniform
  if (type == "unif_norm") {
    x.stdunif <- (x - min(x)) / (max(x) - min(x))
    x.stdunif[x.stdunif == 1] <- 1 - 1e-5
    x.stdunif[x.stdunif == 0] <- 1e-5
    x.stdnorm <- qnorm(x.stdunif)
    return(
      list(x.trans = x.stdnorm)
    )
  } else if (type == "norm_unif") {
    x.stdunif <- pnorm(x)
    x.unif <- x.stdunif * diff(range) + range[1]
    return(
      list(x.trans = x.unif)
    )
  }
}

#### simulation func =============
get.BN.weight <- function(data, bn.model, type = c("ll", "energy"), n_simu = 5000) {
  if (type == "ll") {
    if (ncol(data) > 1) {
      ll_RCT <- as.numeric(logLik(bn.model, data = data))
      ll_sim <- sapply(1:n_simu, function(rr) {
        sim_data <- rbn(bn.model, n = nrow(data))
        ll <- as.numeric(logLik(bn.model, data = sim_data))
      })
    } else {
      ll_RCT <- sum(dnorm(
        data[[names(data)]],
        mean = bn.model[[names(data)]]$coefficients,
        sd = bn.model[[names(data)]]$sd, log = TRUE
      ))
      ll_sim <- sapply(1:n_simu, function(rr) {
        sim_data <- rbn(bn.model, n = nrow(data))
        ll <- sum(dnorm(
          sim_data[[names(data)]],
          mean = bn.model[[names(data)]]$coefficients,
          sd = bn.model[[names(data)]]$sd, log = TRUE
        ))
      })
    }
    tail_prob <- 2 * min(c(mean(ll_sim < ll_RCT), mean(ll_sim > ll_RCT)))
  } else if (type == "energy") {
    to_mat <- function(df) {
      df %>%
        mutate(across(where(is.factor), as.integer)) %>%
        as.matrix()
    }
    X_v <- to_mat(data)
    n_v <- nrow(X_v)
    m_G <- n_simu
    kappa_G <- 0.5

    # Generate m_G samples from BN in the same (transformed) space
    X_tilde <- to_mat(
      data.frame(rbn(bn.model, n = m_G))[, colnames(data), drop = FALSE]
    )

    # Cross-distance matrix (n_v x m_G)
    cross_d <- sqrt(pmax(
      outer(rowSums(X_v^2), rowSums(X_tilde^2), "+") - 2 * tcrossprod(X_v, X_tilde),
      0
    ))

    # Within-sample distance sums; dist() = lower triangle, full sym sum = 2 * sum(dist)
    d_v <- sum(dist(X_v))
    d_tilde <- sum(dist(X_tilde))

    # Empirical energy distance
    E_hat <- (2 / (n_v * m_G)) * sum(cross_d) -
      (2 / n_v^2) * d_v -
      (2 / m_G^2) * d_tilde

    # Pooled scale: mean pairwise Euclidean distance in X_v union X_tilde
    s_hat <- mean(dist(rbind(X_v, X_tilde)))

    # Normalized distance, effective n, compatibility score
    D_hat <- E_hat / s_hat
    n_eff <- n_v * m_G / (n_v + m_G)
    # tail_prob <- 2 * (1 - pnorm(n_eff^kappa_G * D_hat))
    # tail_prob <- min(max(tail_prob, 0), 1) # clamp; can exceed 1 under sampling noise
    # tail_prob <- 1 - E_hat / (E_hat + 2*s_hat)
    tail_prob <- exp(-n_eff^kappa_G * D_hat)
  }

  return(tail_prob)
}

digital.control <- function(rwd.data, exp.all, EHR.data, RCT.data, 
                            synctrl.n, syn.nset, 
                            trt.eff, seed, bn.type = c(1, 2, 3)) {
  # synctrl.n: sample size of synthetic controls in stage 1 (e.g., 100)
  # bn.type = 1 for using RCT + RWD data in structure learning; bn.type = 2 for using RCT data only
  set.seed(seed)
  
  rwd.data <- rwd.data %>% select(-c(treatment))
  exp.data <- RCT.data %>% filter(treatment == 1)
  ctrl.data <- RCT.data %>% filter(treatment == 0)
  EHR.data <- EHR.data %>% select(-c(treatment, S, Y))
  RCT.data.trans <- select(RCT.data, -c(treatment, S, Y))
  
  ### BART prediction for DC: using only RWD
  bart.model1 <- bart2(
    Y ~ ., data = rwd.data, n.samples = 2500, n.chains = 4, keepTrees = TRUE,
    combineChains = T, n.threads = 1, verbose = FALSE 
  )
  
  ### DT of RCT control arm
  y.pred.dt <- predict(
    bart.model1, (ctrl.data %>% select(-c(treatment, S, Y))), type = "ppd"
  ) %>% colMeans()
  
  ### BN structure learning
  if (bn.type == 1) {
    bn.structure.data <- EHR.data
    bn.transdata <- bn.data <- RCT.data %>% select(-c(treatment, S, Y))
  } else if (bn.type == 2) {
    bn.structure.data <- EHR.data
    bn.transdata <- bn.data <- EHR.data
  } else if (bn.type == 3) {
    bn.structure.data <- bind_rows(
      RCT.data %>% select(-c(treatment, S, Y)), rwd.data %>% select(-c(Y))
    )
    bn.transdata <- bn.data <- RCT.data %>% select(-c(treatment, S, Y))
  }
  ### data transformation
  coltype <- sapply(1:ncol(bn.data), function(r) class(data.frame(bn.data)[, r]))
  if (sum(coltype != "numeric") == 3) {
    for (ii in c(1:ncol(bn.data))[coltype == "numeric"]) {
      bn.transdata[, ii] <- data.transform(pull(bn.data[, ii]), type = "unif_norm")
      RCT.data.trans[, ii] <- data.transform(pull(RCT.data.trans[, ii]), type = "unif_norm")
    }
  } else if (sum(coltype != "numeric") == 2) {
    if (any(grepl("X4", colnames(bn.data)))) {
      bn.transdata["X4"] <- data.transform(bn.data$X4, type = "unif_norm")
      RCT.data.trans["X4"] <- data.transform(RCT.data.trans$X4, type = "unif_norm")
    }
    if (any(grepl("X5", colnames(bn.data)))) {
      bn.transdata["X5"] <- data.transform(bn.data$X5, type = "unif_norm") 
      RCT.data.trans["X5"] <- data.transform(RCT.data.trans$X5, type = "unif_norm") 
    }
  }
  
  # bn.struct: learned from current data
  bn.struct <- pc.stable(bn.structure.data %>% data.frame(), undirected = TRUE)
  bn.order <- rev(colnames(bn.transdata))
  if (any(grepl("mis", bn.order))) {
    bn.order <- c(bn.order[-grep("mis", bn.order)], bn.order[grep("mis", bn.order)])
  }
  bn.struct <- pdag2dag(bn.struct, ordering = bn.order)
  
  ### BN parameter learning
  # exp.motbf2: learned from stage 1 exp data & stage 2 RCT
  rwd.data <- rwd.data %>% data.frame()
  if (all(coltype == "numeric")) {
    bn.model <- bn.fit(bn.struct, data = data.frame(bn.transdata), method = "mle-g")
  } else {
    bn.model <- bn.fit(bn.struct, data = data.frame(bn.transdata), method = "mle-cg")
  }
    
  
  ### generate 100 datasets
  # bn.ctrl2
  bn.ctrl2 <- list()
  idx <- 0
  while (idx < syn.nset) {
    set.seed(233 * idx + 1234)
    bn.transtmp <- bn.tmp <- rbn(bn.model, n = synctrl.n) %>% tibble()
    if (sum(coltype != "numeric") == 3) {
      # transform back
      for (ii in c(1:ncol(bn.tmp))[coltype == "numeric"]) {
        bn.transtmp[, ii] <- data.transform(
          pull(bn.tmp[, ii]),
          range = c(min(bn.data[, ii]), max(bn.data[, ii])),
          type = "norm_unif"
        )
      }
      bn.transtmp <- bn.transtmp %>%
        mutate(X5 = factor(X5), X6 = factor(X6), X7 = factor(X7))
    } else if (sum(coltype != "numeric") == 2) {
      if (any(grepl("X4", colnames(bn.data)))) {
        bn.transtmp["X4"] <- data.transform(
          bn.tmp$X4,
          range = c(min(bn.data$X4), max(bn.data$X4)),
          type = "norm_unif"
        )
      }
      if (any(grepl("X5", colnames(bn.data)))) {
        bn.transtmp["X5"] <- data.transform(
          bn.tmp$X5,
          range = c(min(bn.data$X5), max(bn.data$X5)),
          type = "norm_unif"
        )
      }
      bn.transtmp <- bn.transtmp %>% mutate(X6 = factor(X6), X7 = factor(X7))
    }
    bn.ctrl2 <- list.append(bn.ctrl2, bn.transtmp)
    idx <- idx + 1
  }
  
  output <- lapply(1:syn.nset, function(rr) {
    set.seed(10*rr)
    bart.pred <- predict(bart.model1, data.frame(bn.ctrl2[[rr]]), type = "ppd")
    rdm.idx <- sample(1:nrow(bart.pred), ncol(bart.pred), replace = TRUE)
    y.pred.syn1 <- sapply(1:synctrl.n, function(r) bart.pred[rdm.idx[r], r])
    return(list(y.pred.syn1))
  })
  y.pred.syn1 <- sapply(1:syn.nset, function(r) output[[r]][[1]])
  
  return(list(
    # bn.pval = bn.pval,
    RCT.data.trans = RCT.data.trans,
    bn.model = bn.model, 
    DC.groups = bn.ctrl2,
    y.pred.dt = y.pred.dt, 
    y.pred.syn1 = y.pred.syn1
  ))
}
get.mu.pred <- function(y.pred, prior.shrinkage) {
  jagsdata <- list(
    Ngroup = ncol(y.pred),
    ybar_syn = colMeans(y.pred),
    tau_syn = (nrow(y.pred))/apply(y.pred, MARGIN = 2, var),
    tau_pred = median((nrow(y.pred))/apply(y.pred, MARGIN = 2, var))
  )
  jagsmodel <- run.jags(
    model = str_replace(BHM.pred, "prior_to_be_defined", prior.shrinkage),
    monitor = c("ybar_pred"),
    data = jagsdata, n.chains = 4,
    adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
    method = "rjags", plots = FALSE, silent.jags = T,
    inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
  )
  ybar_pred <- as.mcmc.list(jagsmodel, "ybar_pred") %>% as.matrix()
  return(ybar_pred)
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
MAP.func <- function(rawRWD, RCT.data, true.ctrl.s1, exp.all, trueRCT, 
                     var0.ess, prior.shrinkage = "dunif(0, 100)",
                     wt.rho.x, wt.b.x, wt.rho.y, wt.b.y, wt.type, w0.val, 
                     RCT.data.trans, bn.model, DC.groups, y.pred.dt, y.pred.syn1, 
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
  
  bn.pval <- get.BN.weight(data = data.frame(RCT.data.trans), bn.model = bn.model, type = "energy")
  wt.para[2] <- wt.func(x = bn.pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
  
  wt.para <- pmax(wt.para, 1e-12)
  
  if(outcome.type == 1) {
    ### DC_unadj_v2: std
    jagsdata <- list(
      Ngroup = ncol(y.pred.syn1),
      ybar_syn = colMeans(y.pred.syn1),
      tau_syn = (nrow(y.pred.syn1))/apply(y.pred.syn1, MARGIN = 2, var),
      N_RCT = nrow(RCT.data),
      var0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/(nrow(y.pred.syn1)), 
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
      model = str_replace(unadjMAP.normal, "prior_to_be_defined", prior.shrinkage), 
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
      Ngroup = ncol(y.pred.syn1),
      ybar_syn = colMeans(y.pred.syn1),
      tau_syn = (nrow(y.pred.syn1))/apply(y.pred.syn1, MARGIN = 2, var),
      N_RCT = nrow(RCT.data),
      P = sum(grepl("X", colnames(RCT.data))), 
      var0 = max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(RCT.data$Y))/(nrow(y.pred.syn1)), 
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
      model = str_replace(adjMAP.normal, "prior_to_be_defined", prior.shrinkage), 
      monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"), 
      data = jagsdata, n.chains = 4, 
      adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      method = "rjags", plots = FALSE, silent.jags = T, 
      inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    )
    
    res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data, type = outcome.type)
    prob.DCadj1 <- res.tmp["prob"]
    ATE.DCadj1 <- res.tmp["ATE"]
    

    if (methods == "full") {
      ### 1:1 cross-fitting, no replicate
      shuffled_ctrl <- sample(1:nrow(ctrl.data))
      shuffled_dc <- sample(1:nrow(y.pred.syn1))
      shuffled_trt <- sample(1:nrow(exp.all))
      fold_ctrl <- cut(seq_along(shuffled_ctrl), 2, labels = FALSE)
      fold_dc <- cut(seq_along(shuffled_dc), 2, labels = FALSE)
      fold_trt <- cut(seq_along(shuffled_trt), 2, labels = FALSE)

      # get weights from ML models: cross-validation
      wt.para.cf <- matrix(NA, nrow = 2, ncol = 2)
      for (fold in 1:2) {
        valid_idx <- shuffled_ctrl[fold_ctrl != fold]
        idx_rct <- RCT.data %>% mutate(row_id = row_number()) %>%
          group_by(treatment) %>% mutate(ingroup_idx = row_number()) %>%
          filter(
            (treatment == 0 & ingroup_idx %in% shuffled_ctrl[fold_ctrl != fold]) |
              (treatment == 1 & ingroup_idx %in% shuffled_trt[fold_trt != fold])
          ) %>% ungroup() %>% pull(row_id)

        t_stat <- t.test(y.pred.dt[valid_idx], ctrl.data$Y[valid_idx], paired = TRUE)$statistic
        pval <- 2*(1 - pnorm(abs(t_stat) * length(valid_idx)^(-0.01)))
        wt.para.cf[fold, 1] <- wt.func(x = pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
        bn.pval <- get.BN.weight(data = data.frame(RCT.data.trans[idx_rct, ]), bn.model = bn.model, type = "energy")
        wt.para.cf[fold, 2] <- wt.func(x = bn.pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
        wt.para.cf[fold, ] <- pmax(wt.para.cf[fold, ], 1e-12)
      }

      for (fold in 1:2) {
        idx_dc <- shuffled_dc[fold_dc == fold]
        idx_rct <- RCT.data %>% mutate(row_id = row_number()) %>%
          group_by(treatment) %>% mutate(ingroup_idx = row_number()) %>%
          filter(
            (treatment == 0 & ingroup_idx %in% shuffled_ctrl[fold_ctrl == fold]) |
              (treatment == 1 & ingroup_idx %in% shuffled_trt[fold_trt == fold])
          ) %>% ungroup() %>% pull(row_id)

        jagsdata <- list(
          Ngroup = ncol(y.pred.syn1),
          ybar_syn = colMeans(y.pred.syn1[idx_dc, ]),
          tau_syn = (length(idx_dc))/apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var),
          N_RCT = length(idx_rct),
          var0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/(nrow(y.pred.syn1[idx_dc, ])),
          var_dist0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/var0.ess,
          y_rct = RCT.data[idx_rct, ]$Y,
          treatment = RCT.data[idx_rct, ]$treatment,
          w0 = case_when(
            w0.val == -1 ~ min(wt.para.cf[fold, ]),
            w0.val == 2 ~ wt.para.cf[fold, 1],
            w0.val == 3 ~ wt.para.cf[fold, 2],
            TRUE ~ w0.val
          )
        )
        jagsmodel <- run.jags(
          model = str_replace(unadjMAP.normal, "prior_to_be_defined", prior.shrinkage),
          monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"),
          data = jagsdata, n.chains = 4,
          adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
          method = "rjags", plots = FALSE, silent.jags = T,
          inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
        )
        res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data[idx_rct, ], type = outcome.type, cross.fit = TRUE)
        if (fold == 1) {
          ATE.mcmc <- res.tmp$ATE.mcmc
        } else {
          ATE.mcmc <- cbind(ATE.mcmc, res.tmp$ATE.mcmc)
        }
      }
      ATE.mcmc.pool <- rowMeans(ATE.mcmc)
      prob.DCunadj.cf = 2 * min(mean(ATE.mcmc.pool < 0), mean(ATE.mcmc.pool > 0))
      ATE.DCunadj.cf = mean(ATE.mcmc.pool)
      
      ### DC_adj_cf
      # for (fold in 1:2) {
      #   idx_dc <- shuffled_dc[fold_dc == fold]
      #   idx_rct <- RCT.data %>% mutate(row_id = row_number()) %>%
      #     group_by(treatment) %>% mutate(ingroup_idx = row_number()) %>%
      #     filter(
      #       (treatment == 0 & ingroup_idx %in% shuffled_ctrl[fold_ctrl == fold]) |
      #         (treatment == 1 & ingroup_idx %in% shuffled_trt[fold_trt == fold])
      #     ) %>% ungroup() %>% pull(row_id)
      #   jagsdata <- list(
      #     Ngroup = ncol(y.pred.syn1),
      #     ybar_syn = colMeans(y.pred.syn1[idx_dc, ]),
      #     tau_syn = (length(idx_dc))/apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var),
      #     N_RCT = length(idx_rct),
      #     P = sum(grepl("X", colnames(RCT.data))), 
      #     var0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/(nrow(y.pred.syn1[idx_dc, ])), 
      #     var_dist0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/var0.ess,
      #     y_rct = RCT.data[idx_rct, ]$Y, 
      #     treatment = RCT.data[idx_rct, ]$treatment, 
      #     w0 = case_when(
      #       w0.val == -1 ~ min(wt.para.cf[fold, ]),
      #       w0.val == 2 ~ wt.para.cf[fold, 1],
      #       w0.val == 3 ~ wt.para.cf[fold, 2],
      #       TRUE ~ w0.val
      #     ),
      #     X = (apply(as.matrix(RCT.data[idx_rct, ]), c(1, 2), as.numeric)[, grep("X", names(RCT.data)), drop = FALSE])
      #   )
      #   jagsmodel <- run.jags(
      #     model = str_replace(adjMAP.normal, "prior_to_be_defined", prior.shrinkage), 
      #     monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"), 
      #     data = jagsdata, n.chains = 4, 
      #     adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2, 
      #     method = "rjags", plots = FALSE, silent.jags = T,
      #     inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
      #   )
      #   res.tmp <- get.DC.ATE(jagsmodel = jagsmodel, RCT.data = RCT.data[idx_rct, ], type = outcome.type, cross.fit = TRUE)
      #   if (fold == 1) {
      #     ATE.mcmc <- res.tmp$ATE.mcmc
      #   } else {
      #     ATE.mcmc <- cbind(ATE.mcmc, res.tmp$ATE.mcmc)
      #   }
      # }
      # ATE.mcmc.pool <- rowMeans(ATE.mcmc)
      # prob.DCadj.cf = 2 * min(mean(ATE.mcmc.pool < 0), mean(ATE.mcmc.pool > 0))
      # ATE.DCadj.cf = mean(ATE.mcmc.pool)
      
      ### 1:1 cross-fitting, 9 replicate
      wt.para.cf <- matrix(NA, nrow = 2, ncol = 2)
      for (rep in 1:9) {
        shuffled_ctrl <- sample(1:nrow(ctrl.data))
        shuffled_dc <- sample(1:nrow(y.pred.syn1))
        shuffled_trt <- sample(1:nrow(exp.all))
        fold_ctrl <- cut(seq_along(shuffled_ctrl), 2, labels = FALSE)
        fold_dc <- cut(seq_along(shuffled_dc), 2, labels = FALSE)
        fold_trt <- cut(seq_along(shuffled_trt), 2, labels = FALSE)
        for (fold in 1:2) {
          valid_idx <- shuffled_ctrl[fold_ctrl != fold]
          idx_rct <- RCT.data %>% mutate(row_id = row_number()) %>%
            group_by(treatment) %>% mutate(ingroup_idx = row_number()) %>%
            filter(
              (treatment == 0 & ingroup_idx %in% shuffled_ctrl[fold_ctrl != fold]) |
                (treatment == 1 & ingroup_idx %in% shuffled_trt[fold_trt != fold])
            ) %>% ungroup() %>% pull(row_id)

          t_stat <- t.test(y.pred.dt[valid_idx], ctrl.data$Y[valid_idx], paired = TRUE)$statistic
          pval <- 2*(1 - pnorm(abs(t_stat) * length(valid_idx)^(-0.01)))
          wt.para.cf[fold, 1] <- wt.func(x = pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
          bn.pval <- get.BN.weight(data = data.frame(RCT.data.trans[idx_rct, ]), bn.model = bn.model, type = "energy")
          wt.para.cf[fold, 2] <- wt.func(x = bn.pval, ref.stat = (1 - 0.05), rho = wt.rho.y, b = wt.b.y, type = wt.type)
          wt.para.cf[fold, ] <- pmax(wt.para.cf[fold, ], 1e-12)
        }
        for (fold in 1:2) {
          idx_dc <- shuffled_dc[fold_dc == fold]
          idx_rct <- RCT.data %>% mutate(row_id = row_number()) %>%
            group_by(treatment) %>% mutate(ingroup_idx = row_number()) %>%
            filter(
              (treatment == 0 & ingroup_idx %in% shuffled_ctrl[fold_ctrl == fold]) |
                (treatment == 1 & ingroup_idx %in% shuffled_trt[fold_trt == fold])
            ) %>% ungroup() %>% pull(row_id)

          jagsdata_unadj <- list(
            Ngroup = ncol(y.pred.syn1),
            ybar_syn = colMeans(y.pred.syn1[idx_dc, ]),
            tau_syn = (length(idx_dc))/apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var),
            N_RCT = length(idx_rct),
            var0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/(nrow(y.pred.syn1[idx_dc, ])),
            var_dist0 = max(median(apply(y.pred.syn1[idx_dc, ], MARGIN = 2, var)), var(RCT.data[idx_rct, ]$Y))/var0.ess,
            y_rct = RCT.data[idx_rct, ]$Y,
            treatment = RCT.data[idx_rct, ]$treatment,
            w0 = case_when(
              w0.val == -1 ~ min(wt.para.cf[fold, ]),
              w0.val == 2 ~ wt.para.cf[fold, 1],
              w0.val == 3 ~ wt.para.cf[fold, 2],
              TRUE ~ w0.val
            )
          )
          # jagsdata_adj <- c(
          #   jagsdata_unadj,
          #   list(P = sum(grepl("X", colnames(RCT.data))),
          #        X = (apply(as.matrix(RCT.data[idx_rct, ]), c(1, 2), as.numeric)[, grep("X", names(RCT.data)), drop = FALSE]))
          #   )
          jagsmodel_unadj <- run.jags(
            model = str_replace(unadjMAP.normal, "prior_to_be_defined", prior.shrinkage),
            monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"),
            data = jagsdata_unadj, n.chains = 4,
            adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
            method = "rjags", plots = FALSE, silent.jags = T,
            inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
          )
          # jagsmodel_adj <- run.jags(
          #   model = str_replace(adjMAP.normal, "prior_to_be_defined", prior.shrinkage),
          #   monitor = c("mu_ctrl", "beta_trt", "y_rct_pred"),
          #   data = jagsdata_adj, n.chains = 4,
          #   adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
          #   method = "rjags", plots = FALSE, silent.jags = T,
          #   inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
          # )
          res.tmp.unadj <- get.DC.ATE(jagsmodel = jagsmodel_unadj, RCT.data = RCT.data[idx_rct, ], type = outcome.type, cross.fit = TRUE)
          # res.tmp.adj <- get.DC.ATE(jagsmodel = jagsmodel_adj, RCT.data = RCT.data[idx_rct, ], type = outcome.type, cross.fit = TRUE)
          if (fold == 1) {
            ATE.mcmc.unadj <- res.tmp.unadj$ATE.mcmc
            # ATE.mcmc.adj <- res.tmp.adj$ATE.mcmc
          } else {
            ATE.mcmc.unadj <- cbind(ATE.mcmc.unadj, res.tmp.unadj$ATE.mcmc)
            # ATE.mcmc.adj <- cbind(ATE.mcmc.adj, res.tmp.adj$ATE.mcmc)
          }
        }
        if (rep == 1) {
          ATE.pool.unadj <- rowMeans(ATE.mcmc.unadj)
          # ATE.pool.adj <- rowMeans(ATE.mcmc.adj)
        } else {
          ATE.pool.unadj <- cbind(ATE.pool.unadj, rowMeans(ATE.mcmc.unadj))
          # ATE.pool.adj <- cbind(ATE.pool.adj, rowMeans(ATE.mcmc.adj))
        }
      }
      median.idx1 <- which.min(abs(colMeans(ATE.pool.unadj) - median(colMeans(ATE.pool.unadj))))
      prob.DCunadj.cf2 = 2 * min(mean(ATE.pool.unadj[, median.idx1] < 0), mean(ATE.pool.unadj[, median.idx1] > 0))
      ATE.DCunadj.cf2 = mean(ATE.pool.unadj[, median.idx1])
      
      ### 1:1 cross-fitting, 9 replicate end
      
      # median.idx2 <- which.min(abs(colMeans(ATE.pool.adj) - median(colMeans(ATE.pool.adj))))
      # prob.DCadj.cf2 = 2 * min(mean(ATE.pool.adj[, median.idx2] < 0), mean(ATE.pool.adj[, median.idx2] > 0))
      # ATE.DCadj.cf2 = mean(ATE.pool.adj[, median.idx2])
      
    }
    
    
  } else if (outcome.type == 2) {
    jagsdata <- list(
      Ngroup = ncol(y.pred.syn1),
      ysum_syn = colSums(y.pred.syn1),
      n_dc = nrow(y.pred.syn1), 
      N_RCT = nrow(RCT.data),
      var0 = 1/(max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(ctrl.data$Y)) * nrow(y.pred.syn1)), 
      var_dist0 = 1/(max(median(apply(y.pred.syn1, MARGIN = 2, var)), var(ctrl.data$Y)) * var0.ess),
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
      model = str_replace(unadjMAP.binary, "prior_to_be_defined", prior.shrinkage), 
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
    
    prob.DCunadj.cf = NA
    ATE.DCunadj.cf = NA
    
    prob.DCunadj.cf2 = NA
    ATE.DCunadj.cf2 = NA
    
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
        pspp.data = pspp.data, strata.n = ii, borrow.n = nrow(y.pred.syn1), 
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
        MAP_ESS = nrow(y.pred.syn1), S = ii, showESS = FALSE, type = 1, outcome.type = outcome.type
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
    
    ### ATE: PS-SAM (CSD = 0.4); 
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
      ) %>% dplyr::sample_n(nrow(y.pred.syn1)/2, replace = TRUE)
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
    ) %>% dplyr::sample_n(nrow(y.pred.syn1), replace = TRUE)
    
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
    # procova.data <- RCT.data %>% select(-c(S)) %>% mutate(prog.score = prog.score)
    # for (ii in 1:(ncol(procova.data) - 3)) {
    #   if (class(procova.data[[ii]]) == "numeric") {
    #     procova.data[ii] <- procova.data[ii] - mean(pull(procova.data[ii]))
    #   } else if (class(procova.data[[ii]]) == "factor") {
    #     numeric_values <- as.numeric(as.character(procova.data[[ii]]))
    #     procova.data[[ii]] <- numeric_values - mean(numeric_values)
    #   }
    # }
    # procova.data$prog.score <- procova.data$prog.score - mean(procova.data$prog.score)
    # if(outcome.type == 1) {
    #   procova.fit <- lm(Y ~ ., data = procova.data)
    #   ATE.procova <- procova.fit$coefficients["treatment"]
    #   pval.procova <- coef(summary(procova.fit))["treatment", "Pr(>|t|)"]
    # } else if (outcome.type == 2) {
    #   procova.fit <- glm(Y ~ ., data = procova.data, family = binomial)
    #   ATE.procova <- mean(predict(procova.fit, newdata = (procova.data %>% mutate(treatment = 1)), type = "response")) - mean(predict(procova.fit, newdata = (procova.data %>% mutate(treatment = 0)), type = "response"))
    #   pval.procova <- coef(summary(procova.fit))["treatment", "Pr(>|z|)"]
    # }
    
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
      ATE.DCunadj.cf = ATE.DCunadj.cf,
      ATE.DCunadj.cf2 = ATE.DCunadj.cf2,
      ATE.DCadj1 = ATE.DCadj1,
      # ATE.DCtwin = ATE.DCtwin, 
      # ATE.DCtwin2 = ATE.DCtwin2, 
      ATE.anova = ATE.anova,
      ATE.anova2 = ATE.anova2,
      # ATE.anova3 = ATE.anova3,
      ATE.ancova = ATE.ancova
      # ATE.RCTbart = ATE.RCTbart,
      # ATE.pool.lmrobust = ATE.pool.lmrobust,
      # ATE.poolbart = ATE.poolbart,
      # ATE.procova = ATE.procova
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
      prob.DCunadj.cf = prob.DCunadj.cf,
      prob.DCunadj.cf2 = prob.DCunadj.cf2,
      prob.DCadj1 = prob.DCadj1,
      pval.anova = pval.anova,
      pval.anova2 = pval.anova2,
      # pval.anova3 = pval.anova3,
      pval.ancova = pval.ancova
      # pval.RCTbart = pval.RCTbart,
      # pval.pool.lmrobust = pval.pool.lmrobust,
      # pval.poolbart = pval.poolbart,
      # pval.procova = pval.procova
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
      pval = pval
    )
    prob.res <- tibble(
      pval.anova = pval.anova,
      pval.anova2 = pval.anova2,
      pval.ancova = pval.ancova,
      prob.DCunadj2 = prob.DCunadj2,
      prob.DCadj1 = prob.DCadj1
    )
  }
  
  return(list(ATE = ATE.res, Prob = prob.res))
}



### MAIN function ==================================
MAIN.func <- function(rwd.n, exp.n, EHR.n, synctrl.n, trt.eff, bias.c, syn.nset,
                      scenario, var0.ess, prior.shrinkage, 
                      wt.rho.x, wt.b.x, wt.rho.y, wt.b.y, wt.type, w0.val, 
                      sigma.rwdx = 1, sigma.rwd = 1, 
                      sigma.rctx = 1, sigma.rct = 1, rho.rwd = 0.3, 
                      model.type, bias.type, bn.type, outcome.type, seed, rep) {
  # syn.nset: number of generate synthetic datasets
  # var0.ess: sigma0^2 ESS value for the non-informative part in the mixture prior

  ### Data generating
  tmp.data <- prepare.data(
    rwd.n = rwd.n, exp.n = exp.n, EHR.n = EHR.n, 
    trt.eff = trt.eff, bias.c = bias.c,
    syn.nset = syn.nset, scenario = scenario, 
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
    syn.nset = syn.nset, 
    trt.eff = trt.eff, 
    seed = seed, 
    bn.type = bn.type
  )
  ### MAP
  res.s2 <- MAP.func(
    rawRWD = rawRWD, 
    RCT.data = tmp.data$RCT.data,
    true.ctrl.s1 = true.ctrl.s1, 
    exp.all = exp.all, 
    trueRCT = tmp.data$trueRCT, 
    var0.ess = var0.ess,
    prior.shrinkage = prior.shrinkage,
    wt.rho.x = wt.rho.x, 
    wt.b.x = wt.b.x, 
    wt.rho.y = wt.rho.y, 
    wt.b.y = wt.b.y, 
    wt.type = wt.type, 
    w0.val = w0.val, 
    RCT.data.trans = res.s1$RCT.data.trans, 
    bn.model = res.s1$bn.model, 
    DC.groups = res.s1$DC.groups,
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

