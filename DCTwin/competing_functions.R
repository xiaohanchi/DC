ATE.pw <- function(RWD, ctrl.data, exp.all, type = 1) {
  # type: 1 for using trt pts for weighting; 2 for using ctrl pts.
  set.seed(123)
  exp.subset <- exp.all[sample(1:nrow(exp.all), (nrow(exp.all) - nrow(ctrl.data)), replace = FALSE), ]
  
  if (type == 1) {
    ps.data <- bind_rows(
      (RWD %>% mutate(label = 0)),
      (exp.subset %>% mutate(label = 1))
    )
  } else if (type == 2) {
    ps.data <- bind_rows(
      (RWD %>% mutate(label = 0)),
      (ctrl.data %>% mutate(label = 1))
    )
  }
  
  Cname <- names(ps.data)[grep("X", names(ps.data))]
  ## propensity score trimming
  ps.mod <- glm(
    as.formula(paste("label ~ ", paste(Cname, collapse = "+"))),
    data = ps.data, family = binomial
  )
  ps <- predict(ps.mod, newdata = ps.data, type = "response")
  ps.data <- ps.data %>% mutate(ps = ps)
  ps.range <- c(
    min(filter(ps.data, label == 1)$ps),
    max(filter(ps.data, label == 1)$ps)
  )
  ps.data <- ps.data %>%
    filter(ps >= ps.range[1], ps <= ps.range[2]) %>%
    select(-c(ps))
  
  ps.mod <- glm(
    as.formula(paste("label ~ ", paste(Cname, collapse = "+"))),
    data = ps.data, family = binomial
  )
  ps <- predict(ps.mod, newdata = ps.data, type = "response")
  ps.data <- ps.data %>% mutate(wt = ifelse(ps.data$label == 1, 1, ps / (1 - ps)))
  ps.ctrl <- ps.data %>% filter(label == 0)
  
  test.data <- bind_rows(
    exp.all %>% mutate(treatment = 1, wt = 1),
    ps.ctrl %>% select(-c(label)) %>% mutate(treatment = 0, .before = wt),
    ctrl.data %>% mutate(treatment = 0, wt = 1)
  )
  
  pw.fit <- lm(Y ~ treatment, weights = wt, data = test.data)
  
  # Ctrl.mu <- (sum(ps.ctrl$Y * ps.ctrl$wt) + sum(ctrl.data$Y)) / (sum(ps.ctrl$wt) + nrow(ctrl.data))
  # Trt.mu <- mean(exp.all$Y)
  
  ATE <- coef(summary(pw.fit))["treatment", 1]
  pval <- pt(coef(summary(pw.fit))["treatment", 3], pw.fit$df, lower = FALSE)
  return(
    list(ATE = ATE, pval = pval)
  )
}

# ATE.basic <- function(psdata.s1, psdata.full, ctrl.full, n1, N, pi.val = 0.9, L.simu) {
#   basic.ps <- glm(label ~ . - Y, data = psdata.s1, family = binomial()) %>%
#     augment(type.predict = "response", data = psdata.s1) %>%
#     rename("PS" = ".fitted") %>%
#     mutate(logitPS = logit(PS))
#   
#   post.mu <- dist_student_t(
#     df = (nrow(basic.ps) - 1), mu = mean(basic.ps$logitPS), sigma = sqrt(var(basic.ps$logitPS) / nrow(basic.ps))
#   )
#   post.trans.sig2 <- dist_chisq(nrow(basic.ps) - 1)
#   sample.mu <- unlist(generate(post.mu, L.simu))
#   sample.sig2 <- ((nrow(basic.ps) - 1) * var(basic.ps$logitPS)) / unlist(generate(post.trans.sig2, L.simu))
#   
#   pred.Ns <- sapply(1:L.simu, function(r) {
#     set.seed(233 * r)
#     logitPS <- rnorm((N - n1), sample.mu[r], sample.sig2[r])
#     PS <- expit(logitPS)
#     pred.exp <- tibble(label = 1, PS = PS)
#     all.data <- bind_rows(pred.exp, basic.ps[c("label", "PS")])
#     match.pred <- matchit(label ~ PS, all.data,
#                           method = "nearest", distance = all.data$PS,
#                           caliper = 0.2, std.caliper = TRUE, ratio = 1, tol = 1e-10
#     )
#     Ns <- match.data(match.pred) %>%
#       filter(label == 1) %>%
#       nrow()
#     return(Ns)
#   })
#   SynEff <- median(pred.Ns) / N
#   
#   # SC
#   match.it <- matchit(label ~ . - Y,
#                       data = psdata.full,
#                       method = "nearest",
#                       caliper = 0.2, std.caliper = TRUE,
#                       ratio = 1, tol = 1e-10
#   )
#   matchdata <- match.data(match.it)[1:ncol(psdata.full)]
#   
#   if (SynEff >= pi.val) {
#     # SC
#     ATE <- mean(filter(matchdata, label == 1)$Y) - mean(filter(matchdata, label == 0)$Y)
#     pval <- t.test(filter(matchdata, label == 1)$Y, filter(matchdata, label == 0)$Y, alternative = "greater")$p.value
#     sample.size <- N
#   } else {
#     ATE <- mean(filter(matchdata, label == 1)$Y) -
#       mean(c(filter(matchdata, label == 0)$Y, ctrl.full$Y[1:(N - SynEff * N)]))
#     pval <- t.test(
#       filter(matchdata, label == 1)$Y,
#       c(filter(matchdata, label == 0)$Y, ctrl.full$Y[1:(N - SynEff * N)]),
#       alternative = "greater"
#     )$p.value
#     sample.size <- 2 * N - SynEff * N
#   }
#   return(list(ATE = ATE, pval = pval, sample.size = sample.size))
# }


ATE.pspp <- function(pspp.data, strata.n = 5, borrow.n, type = 1, outcome.type) {
  # type = 1 for using all current data; type = 2 for using only current control data
  # outcome.type: "continuous" or "binary"
  if (type == 1) {
    ### Obtain PSs.
    ps <- psrwe_est(
      data.frame(pspp.data),
      v_covs = colnames(pspp.data)[grep("X", colnames(pspp.data))],
      v_grp = "label", cur_grp_level = "1",
      v_arm = "arm", ctl_arm_level = "0",
      ps_method = "logistic", nstrata = strata.n, stra_ctl_only = FALSE
    )
    borrow.res <- psrwe_borrow(ps, total_borrow = borrow.n)
    options(mc.cores = 1)
    .msg <- capture.output({
      suppressWarnings({
        output <- psrwe_powerp(borrow.res,
                               outcome_type = outcome.type,
                               v_outcome    = "Y",
                               seed         = 2333
        )
      })
    })
    ATE <- output$Effect$Overall_Estimate[1, "Mean"]
    prob <- 2 * min(mean(output$Effect$Overall_Samples < 0), mean(output$Effect$Overall_Samples > 0))
    return(
      list(ATE = ATE, Prob = prob)
    )
  } else if (type == 2) {
    ### Obtain PSs.
    ps <- psrwe_est(
      data.frame(pspp.data),
      v_covs = colnames(pspp.data)[grep("X", colnames(pspp.data))],
      v_grp = "label", cur_grp_level = "1",
      ps_method = "logistic", nstrata = strata.n, stra_ctl_only = TRUE
    )
    borrow.res <- psrwe_borrow(ps, total_borrow = borrow.n)
    options(mc.cores = 1)
    .msg <- capture.output({
      suppressWarnings({
        output <- psrwe_powerp(borrow.res,
                               outcome_type = "continuous",
                               v_outcome    = "Y",
                               seed         = 2333
        )
      })
    })
    
    muC <- output$Control$Overall_Estimate[1, "Mean"]
    return(
      list(muC = muC)
    )
  }
}

ATE.pssam <- function(ps.data, ctrl.data, exp.data, sigma = NULL, eff.size, outcome.type){
  ps.data <- data.frame(ps.data)
  prior1 <- PS_prior(
    formula = paste0("label ~ ", paste(names(ps.data)[grep("X", names(ps.data))], collapse = "+")),
    data = ps.data, ps.method = 'Matching', method = 'nearest', 
    outcome = 'Y', study = 'label', treat = 'arm'
    )
  
  wSAM <- SAM_weight(
    if.prior = prior1,
    delta = eff.size,    ## Clinically significant difference
    data = ctrl.data$Y     ## Control arm data
  )
  if(outcome.type == 1) {
    nf.prior <- mixnorm(nf.prior = c(1, 0, sigma), sigma = sigma)
  } else if (outcome.type == 2) {
    nf.prior <- mixbeta(nf.prior = c(1, 1, 1))
  }
  
  SAM.prior <- SAM_prior(if.prior = prior1, nf.prior = nf.prior, weight = wSAM)
  
  post_SAM <- postmix(priormix = SAM.prior, data = ctrl.data$Y)
  post_trt <- postmix(priormix = nf.prior, data = exp.data$Y)
  n_samp <- 20000
  ctrl_samp <- rmix(post_SAM, n_samp)
  trt_samp <- rmix(post_trt, n_samp)

  results <- list(
    Prob = 2 * min(mean((trt_samp - ctrl_samp) < 0), mean((trt_samp - ctrl_samp) > 0)), 
    ATE = mean(trt_samp) - mean(ctrl_samp)
  )
  return(results)
  
}

ATE.pssam.cf <- function(ps.data, ctrl.data.valid, ctrl.data.est, exp.data, sigma = NULL, eff.size, outcome.type){
  ps.data <- data.frame(ps.data)
  prior1 <- PS_prior(
    formula = paste0("label ~ ", paste(names(ps.data)[grep("X", names(ps.data))], collapse = "+")),
    data = ps.data, ps.method = 'Matching', method = 'nearest', 
    outcome = 'Y', study = 'label', treat = 'arm'
  )
  
  wSAM <- SAM_weight(
    if.prior = prior1,
    delta = eff.size,    ## Clinically significant difference
    data = ctrl.data.valid$Y     ## Control arm data
  )
  if(outcome.type == 1) {
    nf.prior <- mixnorm(nf.prior = c(1, 0, sigma), sigma = sigma)
  } else if (outcome.type == 2) {
    nf.prior <- mixbeta(nf.prior = c(1, 1, 1))
  }
  
  SAM.prior <- SAM_prior(if.prior = prior1, nf.prior = nf.prior, weight = wSAM)
  
  post_SAM <- postmix(priormix = SAM.prior, data = ctrl.data.est$Y)
  post_trt <- postmix(priormix = nf.prior, data = exp.data$Y)
  n_samp <- 20000
  ctrl_samp <- rmix(post_SAM, n_samp)
  trt_samp <- rmix(post_trt, n_samp)
  
  results <- list(
    ATE.mcmc = (trt_samp - ctrl_samp)
  )
  return(results)
  
}



### PW-MEM ====================================================
get.part <- function(R, max_cl) {
  ## generate all the possible partitions
  ## and store them in a matrix
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x) {
    length(unique(x)) <= max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  return(part_mat)
}


set.mem.prior <- function(num_study, delta) {
  part <- get.part(R = num_study, max_cl = num_study)
  K <- nrow(part)
  # number of blocks/unique response rate in each partition
  n_bk <- apply(part, 1, function(x) {
    length(unique(x))
  })
  prior <- n_bk^(delta) / sum(n_bk^(delta))
  
  list("part" = part, "prior" = prior)
}

update.part.normal <- function(u, v, prior_part, part) {
  R <- length(u)
  K <- nrow(part)
  
  p <- foreach(k = 1:K, .combine = "cbrob") %do% {
    grp <- unlist(part[k, ])
    U <- aggregate(u, by = list(grp), sum)$x
    V <- aggregate(v, by = list(grp), sum)$x
    # calculate marginal probs m(s_j)
    log_likelood <- as.brob(U^2 / (4 * V) + 1 / 2 * log(pi / V))
    prod(exp(log_likelood)) * prior_part[k]
  }
  post_part <- as.numeric(p / sum(p))
  
  idx <- which.max(post_part)
  
  sim_mat <- part[, 1:(R - 1)] == part[, R]
  
  if (is.null(nrow(sim_mat))) {
    post_sim <- sum(sim_mat * post_part)
  } else {
    post_sim <- colSums(sim_mat * post_part)
  }
  
  return(list(
    "part_hat" = part[idx, ], "phat" = post_part[idx],
    "post_part" = post_part,
    "post_sim" = post_sim
  ))
}


update.part.bin <- function(x, n, prior_part, part, a0 = 1, b0 = 1) {
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K, .combine = "c") %do% {
    grp <- unlist(part[k, ])
    S <- aggregate(x, by = list(grp), sum)$x
    N <- aggregate(n, by = list(grp), sum)$x
    # calculate marginal probs m(s_j)
    prod((beta(a0 + S, b0 + N - S) / beta(a0, b0))) * prior_part[k]
  }
  post_part <- p / sum(p)
  
  idx <- which.max(post_part)
  
  sim_mat <- part[, 1:(R - 1)] == part[, R]
  
  if (is.null(nrow(sim_mat))) {
    post_sim <- sum(sim_mat * post_part)
  } else {
    post_sim <- colSums(sim_mat * post_part)
  }
  
  return(list(
    "part_hat" = part[idx, ], "phat" = post_part[idx],
    "post_part" = post_part,
    "post_sim" = post_sim
  ))
}

bayes.two.normal <- function(tau0, mu0, tau1, mu1,
                             prior0 = c(0.5, 0.5),
                             prior1 = c(0.5, 0.5)) {
  theta0 <- rnorm(1e5, mu0, sqrt(1 / tau0))
  theta1 <- rnorm(1e5, mu1, sqrt(1 / tau1))
  
  prob <- 2 * min(mean(theta1 - theta0 < 0), mean(theta1 - theta0 > 0))
  
  return(prob)
}




bayes.two.prop <- function(y0, n0, y1, n1, 
                           prior0 = c(0.5, 0.5),
                           prior1 = c(0.5, 0.5)) {
  post0 <- c(y0 + prior0[1], n0 - y0 + prior0[2])
  post1 <- c(y1 + prior1[1], n1 - y1 + prior1[2])
  theta0 <- rbeta(1e5, post0[1], post0[2])
  theta1 <- rbeta(1e5, post1[1], post1[2])
  
  prob <- 2 * min(mean(theta1 - theta0 < 0), mean(theta1 - theta0 > 0))
  SD = sd(theta1 - theta0)
  CI <- quantile((theta1 - theta0), probs = c(0.025, 0.975))
  HPD <- HPDinterval(as.mcmc(theta1 - theta0), prob = 0.95)
  return(list(prob = prob, SD = SD, CI = CI, HPD = HPD["var1", ]))
}



ATE.pwmem <- function(RWD, ctrl.data, exp.all, outcome.type) {
  ps.data <- bind_rows(
    (RWD %>% mutate(label = 0)),
    (ctrl.data %>% mutate(label = 1))
  ) %>% data.frame()
  Cname <- names(ps.data)[grep("X", names(ps.data))]
  ## propensity score model
  ps.mod <- glm(
    as.formula(paste("label ~ ", paste(Cname, collapse = "+"))),
    data = ps.data, family = binomial
  )
  ps <- predict(ps.mod, newdata = ps.data, type = "response")
  ps.data$wt <- ifelse(ps.data$label == 1, 1, ps / (1 - ps))
  
  EC <- ps.data[ps.data$label == 0, ]
  IC <- ps.data[ps.data$label == 1, ]
  
  if (outcome.type == 1){
    EC.wybar <- sum(EC$wt * EC$Y) / sum(EC$wt)
    EC.wvar <- sum(EC$wt * (EC$Y - EC.wybar)^2) / sum(EC$wt)
    EC.wu <- sum(EC$wt * EC$Y) / EC.wvar
    EC.wv <- 1 / 2 * sum(EC$wt) / EC.wvar
    
    IC.ybar <- mean(IC$Y)
    IC.var <- mean((IC$Y - IC.ybar)^2)
    IC.u <- sum(IC$Y) / IC.var
    IC.v <- 1 / 2 * nrow(IC) / IC.var
    
    fit0 <- set.mem.prior(num_study = 2, delta = 0)
    mem.fit <- update.part.normal(
      u = c(IC.u, EC.wu),
      v = c(IC.v, EC.wv),
      prior_part = fit0$prior,
      part = fit0$part
    )
    Ctrl.tau <- nrow(IC) / IC.var + (sum(EC$wt) * mem.fit$post_part[1]) / EC.wvar
    Ctrl.mu <- ((nrow(IC) * IC.ybar / IC.var) + (sum(EC$wt) * mem.fit$post_part[1] * EC.wybar / EC.wvar)) / Ctrl.tau
    
    Trt.tau <- nrow(exp.all) / mean((exp.all$Y - mean(exp.all$Y))^2)
    Trt.mu <- (nrow(exp.all) * mean(exp.all$Y) / mean((exp.all$Y - mean(exp.all$Y))^2)) / Trt.tau
    
    ATE <- Trt.mu - Ctrl.mu
    prob <- bayes.two.normal(tau0 = Ctrl.tau, mu0 = Ctrl.mu, tau1 = Trt.tau, mu1 = Trt.mu)
  } else if (outcome.type == 2) {
    EC_wy <- c(EC$wt%*%EC$Y)
    EC_wn <- sum(EC$wt)
    IC_y <- sum(IC$Y)
    IC_n <- nrow(IC)
    
    IT_y <- sum(exp.all$Y)
    IT_n <- nrow(exp.all)
    fit0 <- set.mem.prior(num_study = 2, delta = 0)
    
    mem.fit <- update.part.bin(
      x = c(IC_y, EC_wy),
      n = c(IC_n, EC_wn),
      prior_part = fit0$prior,
      part = fit0$part
    )
    Ctrl.mu <- (0.5 + IC_y + mem.fit$post_part[1] * EC_wy)/(0.5 + 0.5 + IC_n + mem.fit$post_part[1] * EC_wn)
    Trt.mu <- (0.5 + IT_y)/(0.5 + 0.5 + IT_n)
    
    ATE <- Trt.mu - Ctrl.mu
    tmp <- bayes.two.prop(
      y0 = IC_y + EC_wy*mem.fit$post_part[1], 
      n0 = IC_n + EC_wn*mem.fit$post_part[1], 
      y1 = IT_y, n1 = IT_n
    )
    prob <- tmp$prob
  }
  
  return(
    list(ATE = ATE, Prob = prob)
  )
}


### PS-MAP ====================================================

MAP_Prior <- "
model {  
 ## Likelihoods ##
	for(i in 1:n.hist) {
		ybar[i] ~ dnorm(theta[i], tau.hat[i])
		theta[i] ~ dnorm(mu,tau)
	}

	precision_heterogeneity <- 1/(std_heterogeneity*std_heterogeneity)

	#Half-normal prior for std of study effects
	sigma_heterogeneity ~ dnorm(.00001, precision_heterogeneity)

	#Precision of study effects
	tau <- 1/(sigma_heterogeneity*sigma_heterogeneity)

 
	mu ~ dnorm(0, 0.0001)

	theta.new ~ dnorm(mu, tau)


}

"

MAP_Prior_bin <- "
model {
  ## Likelihoods ##
  for (i in 1:N) {
    logit(p[i]) <- delta[i]
    delta[i] ~ dnorm(beta0, tau)
    resp[i] ~ dbin(p[i],n[i])
  }

  precision_heterogeneity <- 1/(std_heterogeneity*std_heterogeneity)

  #Half-normal prior for std of study effects
  sigma_heterogeneity ~ dnorm(.00001, precision_heterogeneity)
  #Precision of study effects
  tau <- 1/(sigma_heterogeneity*sigma_heterogeneity)

  beta0 ~ dnorm(0, 0.0001)



  # To get ESS
  logit_p.pred ~  dnorm(beta0, tau)

  p.pred <- 1 / (1 + exp(-logit_p.pred))

}

"


rwePS <- function(data, ps.fml = NULL,
                  v.grp   = "group",
                  v.arm   = "arm",
                  v.covs  = "V1",
                  d1.arm  = NULL,
                  d1.grp  = 1,
                  nstrata = 5, ...) {
  
  dnames <- colnames(data);
  stopifnot(v.grp %in% dnames);
  
  ## generate formula
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(v.grp, "~",
                               paste(v.covs, collapse = "+"),
                               sep = ""))
  
  ## d1 index will be kept in the results
  d1.inx   <- d1.grp == data[[v.grp]];
  keep.inx <- which(d1.inx);
  
  ## for 2-arm studies only
  if (!is.null(d1.arm))
    d1.inx <- d1.inx & d1.arm == data[[v.arm]];
  
  ## get ps
  all.ps  <- get_ps(data, ps.fml = ps.fml, ...);
  D1.ps   <- all.ps[which(d1.inx)];
  
  
  ## add columns to data
  grp     <- rep(1, nrow(data));
  grp[which(data[[v.grp]] != d1.grp)] <- 0;
  
  data[["_ps_"]]     <- all.ps;
  data[["_grp_"]]    <- grp;
  data[["_arm_"]]    <- data[[v.arm]];
  
  
  ## stratification
  if (nstrata > 0) {
    strata  <- rweCut(D1.ps, all.ps, breaks = nstrata, keep.inx = keep.inx);
    data[["_strata_"]] <- strata;
  }
  
  ## return
  rst <- list(data    = data,
              ps.fml  = ps.fml,
              nstrata = nstrata);
  class(rst) <- get.rwe.class("DWITHPS");
  
  rst
}


get.rwe.class <- function(c.str = c("DWITHPS", "PSDIST", "D_GPS", "GPSDIST")) {
  c.str <- match.arg(c.str);
  switch(c.str,
         DWITHPS  = "RWE_DWITHPS",
         PSDIST   = "RWE_PSDIST",
         D_GPS    = "RWE_D_GPS",
         GPSDIST  = "RWE_GPSDIST")
}

## compute propensity scores
get_ps <- function(dta, ps.fml, type = c("logistic", "randomforest"),
                   ntree = 5000,
                   ..., grp = NULL, ps.cov = NULL) {
  
  type <- match.arg(type);
  
  ## generate formula
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"),
                               sep=""));
  
  ## identify grp if passed from formula
  grp <- all.vars(ps.fml)[1];
  
  ## fit model
  switch(type,
         logistic = {
           glm.fit <- glm(ps.fml, family=binomial, data=dta, ...);
           est.ps <- glm.fit$fitted;
         },
         randomforest = {
           dta[[grp]] <- as.factor(dta[[grp]]);
           rf.fit     <- randomForest(ps.fml, data = dta,
                                      ntree = ntree, ...);
           est.ps     <- predict(rf.fit, type = "prob")[,2];
         });
  est.ps
}

rweCut <- function(x, y = x, breaks = 5, keep.inx = NULL) {
  cuts    <- quantile(x, seq(0, 1,length = breaks+1));
  cuts[1] <- cuts[1] - 0.001;
  rst     <- rep(NA, length(y));
  for (i in 2:length(cuts)) {
    inx      <- which(y > cuts[i-1] & y <= cuts[i]);
    rst[inx] <- i-1;
  }
  
  if (!is.null(keep.inx)) {
    inx <- which(y[keep.inx] <= cuts[1]);
    if (0 < length(inx)) {
      rst[keep.inx[inx]] <- 1;
    }
    
    inx <- which(y[keep.inx] > cuts[length(cuts)]);
    if (0 < length(inx)) {
      rst[keep.inx[inx]] <- length(cuts) - 1;
    }
  }
  
  rst
}

rwePSDist <- function(data.withps, n.bins = 10, min.n0 = 10,
                      type = c("ovl", "kl"), d1.arm = NULL, ...) {
  
  f.narm <- function(inx, dataps) {
    
    if (is.null(dataps[["_arm_"]]))
      return(c(length(inx), 0,0));
    
    n0 <- length(which(0 == dataps[inx, "_arm_"]));
    n1 <- length(which(1 == dataps[inx, "_arm_"]));
    
    c(length(inx), n0 ,n1);
  }
  
  stopifnot(inherits(data.withps,
                     what = get.rwe.class("DWITHPS")));
  
  type <- match.arg(type);
  
  dataps   <- data.withps$data;
  nstrata  <- data.withps$nstrata;
  rst      <- NULL;
  for (i in 1:nstrata) {
    
    inx.ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]];
    inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
    n0.01   <- f.narm(which(inx.ps0), dataps);
    n1.01   <- f.narm(which(inx.ps1), dataps);
    
    if (!is.null(d1.arm) & !is.null(dataps[["_arm_"]])) {
      inx.ps0 <- inx.ps0 & d1.arm == dataps[["_arm_"]];
      inx.ps1 <- inx.ps1 & d1.arm == dataps[["_arm_"]];
    }
    
    ps0 <- dataps[which(inx.ps0), "_ps_"];
    ps1 <- dataps[which(inx.ps1), "_ps_"];
    
    if (0 == length(ps0) | 0 == length(ps1))
      warning("No samples in strata");
    
    if (any(is.na(c(ps0, ps1))))
      warning("NA found in propensity scores in a strata");
    
    if (length(ps0) < min.n0) {
      warning("Not enough data in the external data in the current stratum.
                     External data ignored.");
      cur.dist <- 0;
    } else {
      cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...);
    }
    
    rst <- rbind(rst, c(i, n0.01, n1.01, cur.dist));
  }
  
  ## overall
  inx.tot.ps0 <- which(0 == dataps[["_grp_"]]);
  inx.tot.ps1 <- which(1 == dataps[["_grp_"]]);
  n0.tot.01   <- f.narm(inx.tot.ps0, dataps);
  n1.tot.01   <- f.narm(inx.tot.ps1, dataps);
  
  ps0         <- dataps[inx.tot.ps0, "_ps_"];
  ps1         <- dataps[inx.tot.ps1, "_ps_"];
  all.dist    <- rweDist(ps0, ps1, n.bins = nstrata*n.bins, type = type, ...);
  rst         <- rbind(rst, c(0, n0.tot.01, n1.tot.01, all.dist));
  
  
  colnames(rst) <- c("Strata", "N0", "N00", "N01", "N1", "N10", "N11", "Dist");
  rst           <- data.frame(rst);
  class(rst)    <- append(get.rwe.class("PSDIST"), class(rst));
  
  rst
}

rweDist <- function(sample.F0, sample.F1, n.bins = 10, type = c("ovl", "kl"), epsilon = 10^-6) {
  
  type     <- match.arg(type);
  
  smps     <- c(sample.F0, sample.F1);
  n0       <- length(sample.F0);
  n1       <- length(sample.F1);
  
  if (0 == n0 | 0 == n1)
    return(c(n0, n1, NA));
  
  if (1 == length(unique(smps))) {
    cut.smps <- rep(1, n0+n1)
    n.bins   <- 1;
    warning("Distributions for computing distances are degenerate.",
            call. = FALSE);
  } else {
    cut.smps <- rweCut(smps, breaks = n.bins);
  }
  
  rst <- 0;
  for (j in 1:n.bins) {
    n0.j <- length(which(j == cut.smps[1:n0]));
    n1.j <- length(which(j == cut.smps[(n0+1):(n0+n1)]));
    
    rst  <- rst + switch(type,
                         kl = {ep0  <- (n0.j+epsilon)/(n0 + epsilon * n.bins);
                         ep1  <- (n1.j+epsilon)/(n1 + epsilon * n.bins);
                         ep1 * log(ep1/ep0)},
                         ovl = min(n0.j/n0, n1.j/n1));
  }
  
  if ("kl" == type)
    rst <- 1/(1+rst);
  
  rst;
}

PS_MAP.fit <- function(tau.init, target.ESS, sigma, n.cur, ybar.hist, SE.hist, overlap_coefficients, niter, lim, data.indiv) {
  # Obtain prior samples when using the PS-MAP prior with the specified target ESS
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size
  # sigma: the fixed reference scale in ess() for approximating mixture of normal
  # n.cur: number of subjects in each stratum of the current trial
  # ybar.hist: mean in each PS stratum from the historical trials
  # SE.hist: std error in each PS stratum from the historical trials
  # overlap_coefficients: overlap coefficients in each PS stratum
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  
  
  n.strata <- length(ybar.hist)
  total.cur <- sum(n.cur)
  WS1 <- n.cur / total.cur
  ess.res <- c(0)
  tau.res <- c(0)
  direction <- 0
  tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.init # hyperparameter for HN(.)
  stop <- FALSE
  low <- lim[1]
  high <- lim[2]
  tau.scale <- tau.init
  
  while (stop == FALSE) {
    theta.pred <- list()
    # ess.res.str=rep(0,n.strata)
    
    # stop the loop if the ESS does not converge but the lower and upper bound are very close
    # which means too much variation in approximation using conjugate priors with small change in the hyper-parameter
    if ((high - low <= 0.01)) {
      stop <- TRUE
      tau.scale <- tau.res[(which.min(abs(ess.res - target.ESS)[-1]) + 1)]
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    }
    
    # estimate stratum-specific MAP prior
    for (i in 1:n.strata) {
      dataTemp <- list(
        "ybar" = ybar.hist[i],
        "n.hist" = 1,
        "std_heterogeneity" = tau.prior[i],
        "tau.hat" = (1 / SE.hist[i]^2)
      )
      
      model <-
        jags.model(
          file = textConnection(MAP_Prior),
          data = dataTemp,
          n.chains = 1,
          n.adapt = 0.2 * niter,
          quiet = TRUE
        )
      update(model, n.iter = 0.2 * niter, progress.bar = "none") # burn in
      MAP_model <-
        coda.samples(model,
                     variable.names = "theta.new",
                     thin = 5,
                     n.iter = niter, progress.bar = "none"
        )
      
      theta.pred[[i]] <- c(MAP_model[[1]][, 1])
      
      # mix.res <- automixfit(theta.pred[[i]], Nc=1:4, thresh=0, type="norm")
      # ess.res.str[i]=ess( mix.res , method="elir",sigma=sigma)
    }
    
    
    theta.pred <- do.call(cbind, theta.pred)
    theta <- WS1 %*% t(theta.pred) # overall prior as a weighted average of stratum-specific prior
    
    mix.res <- automixfit(theta[1, ], Nc = 1:4, thresh = 0, type = "norm") # approximation as mixture of normal
    ess.res <- c(ess.res, ess(mix.res, method = "elir", sigma = sigma)) # calculate the ESS
    tau.res <- c(tau.res, tau.scale)
    
    
    # Binary search for the value of tau.scale
    if (abs(ess.res[length(ess.res)] - target.ESS) < 5) { # 5 is the threshold value
      direction <- 0
      stop <- TRUE
    } else if (ess.res[length(ess.res)] > target.ESS) {
      direction <- -1
      low <- tau.scale
      tau.scale <- (low + high) / 2
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    } else {
      direction <- 1
      high <- tau.scale
      tau.scale <- (low + high) / 2
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    }
  }
  
  posterior.indiv <- postmix(mix.res, data.indiv) # calculate the posterior given prior and current data
  results <- list(mix.res, ess.res[length(ess.res)], posterior.indiv, tau.res[length(tau.res)])
  return(results)
}



PS_MAP.fit.bin <- function(tau.init, target.ESS, n.cur, ybar.hist, n.hist, overlap_coefficients, niter, lim, data.indiv) {
  # Obtain prior samples when using the PS-MAP prior with the specified target ESS
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size
  
  # n.cur: number of subjects in each stratum from the current trial
  # ybar.hist: mean in each PS stratum from the historical trials
  # n.hist: number of subjects in each stratum from the historical trials
  # overlap_coefficients: overlap coefficients in each PS stratum
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  
  
  n.strata <- length(ybar.hist)
  total.cur <- sum(n.cur)
  WS1 <- n.cur / total.cur
  ess.res <- c(0)
  tau.res <- c(0)
  direction <- 0
  
  tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.init # hyperparameter for HN(.)
  stop <- FALSE
  low <- lim[1]
  high <- lim[2]
  tau.scale <- tau.init
  while (stop == FALSE) {
    theta.pred <- list()
    # ess.res.str=rep(0,n.strata)
    
    # stop the loop if the ESS does not converge but the lower and upper bound are very close
    # which means too much variation in approximation using conjugate priors with small change in the hyper-parameter
    if ((high - low <= 0.01)) {
      stop <- TRUE
      tau.scale <- tau.res[which.min(abs(ess.res - target.ESS))]
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    }
    
    
    # estimate stratum-specific MAP prior
    for (i in 1:n.strata) {
      dataTemp <- list(
        "resp" = ybar.hist[i],
        "n" = n.hist[i],
        N = 1
      )
      dataTemp$std_heterogeneity <- tau.prior[i]
      
      
      model <-
        jags.model(
          file = textConnection(MAP_Prior_bin),
          data = dataTemp,
          n.chains = 1,
          n.adapt = 0.2 * niter,
          quiet = TRUE
        )
      update(model, n.iter = 0.2 * niter, progress.bar = "none") # burn in
      MAP_model <-
        coda.samples(model,
                     variable.names = "p.pred",
                     thin = 5,
                     n.iter = niter, progress.bar = "none"
        )
      
      theta.pred[[i]] <- c(MAP_model[[1]][, 1])
      
      # mix.res <- automixfit(theta.pred[[i]], Nc=1:4, thresh=0, type="beta")
      # ess.res.str[i]=ess( mix.res , method="elir")
    }
    
    
    theta.pred <- do.call(cbind, theta.pred)
    theta <- WS1 %*% t(theta.pred) # overall prior as a weighted average of stratum-specific prior
    
    mix.res <- automixfit(theta[1, ], Nc = 1:5, thresh = 0, type = "beta") # approximation as mixture of beta
    
    ess.res <- c(ess.res, ess(mix.res, method = "elir")) # calculate ESS
    tau.res <- c(tau.res, tau.scale)
    
    # Binary search for the value of tau.scale
    if (abs(ess.res[length(ess.res)] - target.ESS) < 5) {
      direction <- 0
      stop <- TRUE
    } else if (ess.res[length(ess.res)] > target.ESS) {
      direction <- -1
      low <- tau.scale
      tau.scale <- (low + high) / 2
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    } else {
      direction <- 1
      high <- tau.scale
      tau.scale <- (low + high) / 2
      tau.prior <- (min(overlap_coefficients) / overlap_coefficients) * tau.scale
    }
  }
  
  
  # calculate the posterior given prior and current data
  posterior.indiv <- postmix(mix.res, n = length(data.indiv), r = sum(data.indiv))
  results <- list(mix.res, ess.res[length(ess.res)], posterior.indiv, tau.res[length(tau.res)])
  return(results)
}



PS_MAP.fit2 <- function(ybar.trt, sigma.hat.trt, n.trt, tau.init, target.ESS, sigma, n.cur, ybar.hist, SE.hist, overlap_coefficients,
                        lim, data.indiv, prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1, niter) {
  # Fits the PS-MAP prior with target ESS and returns a posterior MCMC object of the treatment effect
  
  # ybar.trt/sigma.hat.trt/n.trt: the sample mean, MLE variance (over n not n-1) and sample size for treatment
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size
  # sigma: the fixed reference scale in ess() for approximating mixture of normal
  # n.cur: number of subjects in each stratum of the current trial
  # ybar.hist: mean from each of the PS strata
  # SE.hist: std error from each of the PS strata
  # overlap_coefficients: overlap coefficients in each PS stratum
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  # prior.trt.var: hyperparameter for the normal prior on treatment mean
  # alpha.sigma/beta.sigma: parameters for the IG prior on treatment/control variance
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  
  
  
  
  burn.in <- (0.2) * niter
  # Create vectors to store parameters
  mu.trt <- rep(NA, niter + burn.in)
  var.trt <- rep(NA, niter + burn.in)
  mu <- rep(NA, niter + burn.in)
  var.control <- rep(NA, niter + burn.in)
  diff.effect <- rep(NA, niter + burn.in)
  
  # Specify hyper-parameters
  prior.var <- prior.trt.var
  alpha.var <- alpha.sigma
  beta.var <- beta.sigma
  
  
  
  # Intialize the parameters
  mu.trt[1] <- rnorm(1, sd = sqrt(prior.var))
  var.trt[1] <- rinvgamma(n = 1, shape = alpha.var, rate = beta.var)
  mu[1] <- rnorm(1, sd = sqrt(prior.var)) # Hyper-prior for control mean is actually jefferys prior i.e propto constant
  diff.effect[1] <- mu.trt[1] - mu[1]
  
  for (i in 2:(niter + burn.in)) {
    # Update treatment mean
    mean.mu.trt <- (n.trt * prior.var * ybar.trt) / (n.trt * prior.var + var.trt[i - 1])
    var.mu.trt <- (var.trt[i - 1] * prior.var) / (n.trt * prior.var + var.trt[i - 1])
    mu.trt[i] <- rnorm(n = 1, mean = mean.mu.trt, sd = sqrt(var.mu.trt))
    
    # Update treatment variance
    a.trt <- n.trt / 2 + alpha.var
    b.trt <- (n.trt / 2) * (sigma.hat.trt + ((ybar.trt - mu.trt[i])^2)) + beta.var
    var.trt[i] <- rinvgamma(n = 1, shape = a.trt, rate = b.trt)
  }
  
  
  # Obtain the posterior for the control group
  PS_MAP_fit <- PS_MAP.fit(
    tau.init, target.ESS, sigma, n.cur, ybar.hist, SE.hist, overlap_coefficients,
    niter = 10000, lim, data.indiv
  )
  
  PS_MAP_ESS <- PS_MAP_fit[[2]]
  PS_MAP_post <- PS_MAP_fit[[3]]
  tau_scale_limits <- PS_MAP_fit[[4]]
  
  mu.con <- rmix(PS_MAP_post, n = niter)
  mu.con <- mu.con[seq(1, length(mu.con), by = 5)]
  mu.trt <- mu.trt[-(1:burn.in)]
  mu.trt <- mu.trt[seq(1, length(mu.trt), by = 5)]
  diff.effect <- mu.trt - mu.con[1:length(mu.trt)]
  
  # PS_MAP.model = mcmc(data.frame(mu = mu.con[1:length(mu.trt)],mu.trt = mu.trt,effect = diff.effect))
  PS_MAP.model <- mcmc(data.frame(effect = diff.effect))
  
  fit.result <- list(mcmc.list(PS_MAP.model), PS_MAP_ESS, tau_scale_limits)
  return(fit.result)
}

PS_MAP.fit2.bin <- function(ybar.trt, n.trt, tau.init, target.ESS, n.cur, ybar.hist, n.hist, overlap_coefficients, lim, data.indiv, alpha0 = 1, beta0 = 1, niter) {
  # Fits the PS-MAP prior with target ESS and returns a posterior MCMC object of the treatment effect
  # ybar.trt/n.trt: the sample mean and sample size for treatment
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size
  # n.cur: number of subjects in each stratum from the current trial
  # ybar.hist: mean in each PS stratum from the historical trials
  # n.hist: number of subjects in each stratum from the historical trials
  # overlap_coefficients: overlap coefficients in each PS stratum
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  # alpha0/beta0: parameters for the beta prior on treatment mean
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  
  
  burn.in <- (0.2) * niter
  # Create vectors to store parameters
  mu.trt <- rbeta(niter + burn.in, ybar.trt + alpha0, n.trt - ybar.trt + beta0)
  
  PS_MAP_fit <- PS_MAP.fit.bin(
    tau.init, target.ESS, n.cur, ybar.hist, n.hist, overlap_coefficients,
    niter, lim, data.indiv
  )
  
  PS_MAP_ESS <- PS_MAP_fit[[2]]
  PS_MAP_post <- PS_MAP_fit[[3]]
  tau_scale_limits <- PS_MAP_fit[[4]]
  mu.con <- rmix(PS_MAP_post, n = niter)
  mu.con <- mu.con[seq(1, length(mu.con), by = 5)]
  mu.trt <- mu.trt[-(1:burn.in)]
  mu.trt <- mu.trt[seq(1, length(mu.trt), by = 5)]
  diff.effect <- mu.trt - mu.con[1:length(mu.trt)]
  
  PS_MAP.model <- mcmc(data.frame(effect = diff.effect))
  
  fit.result <- list(mcmc.list(PS_MAP.model), PS_MAP_ESS, tau_scale_limits)
  return(fit.result)
}



PS_regroup <- function(X.hist.combine, Y.hist.combine, X.cur.con, Y.cur.con, X.cur.trt, p, S = 5, type = 1) {
  # Estimate propensity scores and perform stratification, returns stratified data and stratum summary statistics
  # X.hist.combine: matrix of patient covariates, each row is a subject, historical studies
  # Y.hist.combine: vector of patient outcome,historical studies
  # X.cur.con: matrix of patient covariates, each row is a subject, current control arm
  # Y.cur.con: vector of patient outcome, current control arm
  # p: dimension of covariates
  # S : number of strata
  # type: 
  
  if (type == 1) {
    Y.combine <- rbind(Y.hist.combine, Y.cur.con)
    dta <- rbind(X.hist.combine, X.cur.con, X.cur.trt)
    covs <- NULL
    for (k in 1:p) {
      covs[k] <- paste("V", k, sep = "")
    }
    colnames(dta) <- covs
    group <- c(
      rep(0, dim(X.hist.combine)[1]), 
      rep(1, (dim(X.cur.con)[1] + dim(X.cur.trt)[1]))
    ) # historical data =0; current control=1
    arm <- c(
      rep(0, dim(X.hist.combine)[1]), 
      rep(0, dim(X.cur.con)[1]), 
      rep(1, dim(X.cur.trt)[1])
    )
    dta <- cbind(group, arm, dta)
    dta <- as.data.frame(dta)
    ## estimate PS by logistic regression
    ana.ps <-
      rwePS(
        dta,
        v.grp = "group",
        v.arm   = "arm",
        v.covs = covs,
        d1.arm = 0,
        d1.grp  = 1,
        nstrata = S,
        type = "logistic"
      )
    ## get overlapping coefficients
    ana.ovl <- rwePSDist(ana.ps, d1.arm = 0)
    
  } else if (type == 2){
    Y.combine <- rbind(Y.hist.combine, Y.cur.con)
    dta <- rbind(X.hist.combine, X.cur.con)
    covs <- NULL
    for (k in 1:p) {
      covs[k] <- paste("V", k, sep = "")
    }
    colnames(dta) <- covs
    group <- c(
      rep(0, dim(X.hist.combine)[1]), 
      rep(1, dim(X.cur.con)[1])
    ) # historical data =0; current control=1

    dta <- cbind(group, dta)
    dta <- as.data.frame(dta)
    ## estimate PS by logistic regression
    ana.ps <-
      rwePS(
        dta,
        v.grp = "group",
        v.covs = covs,
        nstrata = S,
        type = "logistic"
      )
    ## get overlapping coefficients
    ana.ovl <- rwePSDist(ana.ps)
  }
  
  ## get rs
  rS <- ana.ovl$Dist[1:S]
  
  X.hist.stratum <- list()
  y.hist.stratum <- list()
  n.hist.stratum <- NULL
  X.cur.stratum <- list()
  y.cur.stratum <- list()
  n.cur.stratum <- NULL
  for (i in 1:S) {
    index <- which((ana.ps$data$group == 0) & (ana.ps$data$"_strata_" == i))
    n.hist.stratum[i] <- length(index)
    # X.hist.stratum[[i]] <- ana.ps$data[index, 2:11]
    y.hist.stratum[[i]] <- Y.combine[index]
    
    indey2 <- which((ana.ps$data$group == 1) & (ana.ps$data$"_strata_" == i))
    n.cur.stratum[i] <- length(indey2)
    # X.cur.stratum[[i]] <- ana.ps$data[indey2, 2:11]
    y.cur.stratum[[i]] <- Y.combine[indey2]
  }
  
  stratum.ybar <- sapply(y.hist.stratum, FUN = mean)
  SE.stratum <- sapply(y.hist.stratum, FUN = sd) / sqrt(n.hist.stratum)
  # stratum.ybar.cur <- sapply(y.cur.stratum, FUN = mean)
  # SE.stratum.cur <- sapply(y.cur.stratum, FUN = sd) / sqrt(n.cur.stratum)
  # sigma.stratum.cur <- sapply(y.cur.stratum, FUN = sd) * (n.cur.stratum - 1) / n.cur.stratum
  out <- list(
    rS = rS, X.hist = X.hist.stratum, Y.hist = y.hist.stratum, Ybar.hist = stratum.ybar,
    SE.hist = SE.stratum, 
    # Ybar.cur = stratum.ybar.cur, SE.cur = SE.stratum.cur, sigma.cur = sigma.stratum.cur,
    n.hist = n.hist.stratum, n.cur = n.cur.stratum, PS = ana.ps
  )
  return(out)
}

PS_regroup_bin <- function(X.hist.combine, Y.hist.combine, X.cur.con, Y.cur.con, X.cur.trt, p, S = 5, type = 1) {
  # Estimate propensity scores and perform stratification, returns stratified data and stratum summary statistics
  # X.hist.combine: matrix of patient covariates, each row is a subject, historical studies
  # Y.hist.combine: vector of patient outcome,historical studies
  # X.cur.con: matrix of patient covariates, each row is a subject, current control arm
  # Y.cur.con: vector of patient outcome, current control arm
  # p: dimension of covariates
  
  if (type == 1) {
    Y.combine <- rbind(Y.hist.combine, Y.cur.con)
    dta <- rbind(X.hist.combine, X.cur.con, X.cur.trt)
    covs <- NULL
    for (k in 1:p) {
      covs[k] <- paste("V", k, sep = "")
    }
    colnames(dta) <- covs
    group <- c(
      rep(0, dim(X.hist.combine)[1]),
      rep(1, (dim(X.cur.con)[1] + dim(X.cur.trt)[1]))
    ) # historical data =0; current control=1
    arm <- c(
      rep(0, dim(X.hist.combine)[1]),
      rep(0, dim(X.cur.con)[1]),
      rep(1, dim(X.cur.trt)[1])
    )
    dta <- cbind(group, arm, dta)
    dta <- as.data.frame(dta)
    ## estimate PS by logistic regression
    ana.ps <- rwePS(
      dta,
      v.grp = "group",
      v.arm = "arm",
      v.covs = covs,
      d1.arm = 0,
      d1.grp = 1,
      nstrata = S,
      type = "logistic"
    )
    
    ## get overlapping coefficients
    ana.ovl <- rwePSDist(ana.ps, d1.arm = 0)
  } else if (type == 2) {
    Y.combine <- rbind(Y.hist.combine, Y.cur.con)
    dta <- rbind(X.hist.combine, X.cur.con)
    covs <- NULL
    for (k in 1:p) {
      covs[k] <- paste("V", k, sep = "")
    }
    colnames(dta) <- covs
    group <- c(
      rep(0, dim(X.hist.combine)[1]),
      rep(1, dim(X.cur.con)[1])
    ) # historical data =0; current control=1
    dta <- cbind(group, dta)
    dta <- as.data.frame(dta)
    ## estimate PS by logistic regression
    ana.ps <-
      rwePS(
        dta,
        v.grp = "group",
        v.covs = covs,
        nstrata = S,
        type = "logistic"
      )
    
    ## get overlapping coefficients
    ana.ovl <- rwePSDist(ana.ps)
  }
  
  
  
  ## get rs
  rS <- ana.ovl$Dist[1:S]
  
  X.hist.stratum <- list()
  y.hist.stratum <- list()
  n.hist.stratum <- NULL
  X.cur.stratum <- list()
  y.cur.stratum <- list()
  n.cur.stratum <- NULL
  for (i in 1:S) {
    index <- which((ana.ps$data$group == 0) & (ana.ps$data$"_strata_" == i))
    n.hist.stratum[i] <- length(index)
    # X.hist.stratum [[i]] <- ana.ps$data[index,2:(p+1)]
    y.hist.stratum[[i]] <- Y.combine[index]
    
    indey2 <- which((ana.ps$data$group == 1) & (ana.ps$data$"_strata_" == i) & (ana.ps$data$arm == 0))
    n.cur.stratum[i] <- length(indey2)
    # X.cur.stratum [[i]] <- ana.ps$data[indey2,2:(p+1)]
    y.cur.stratum[[i]] <- Y.combine[indey2]
  }
  
  stratum.ybar <- sapply(y.hist.stratum, FUN = sum)
  
  stratum.ybar.cur <- sapply(y.cur.stratum, FUN = sum)
  
  
  out <- list(
    rS = rS, X.hist = X.hist.stratum, Y.hist = y.hist.stratum,
    Ybar.hist = stratum.ybar, Ybar.cur = stratum.ybar.cur,
    n.hist = n.hist.stratum, n.cur = n.cur.stratum, PS = ana.ps
  )
  return(out)
}


run.psmap <- function(hist.data, current.ctrl, current.trt, MAP_ESS, S, outcome.type, 
                      limits = c(0.001, 2), getlimits = FALSE, type = 1, niter = 1e5){
  
  # type = 1 for using all current data in PS; type = 2 for using control data in PS
  
  p <- grep("X", colnames(hist.data)) %>% length()
  n.trt <- nrow(current.trt)
  
  X.hist.combine <- hist.data %>% select(-c(Y)) %>% data.frame() %>% as.matrix()
  Y.hist.combine <- hist.data %>% select(c(Y)) %>% data.frame() %>% as.matrix()
  X.cur.con <- current.ctrl %>% select(-c(Y)) %>% data.frame() %>% as.matrix()
  y.control <- current.ctrl %>% select(c(Y)) %>% data.frame() %>% as.matrix()
  X.cur.trt <- current.trt %>% select(-c(Y)) %>% data.frame() %>% as.matrix()
  mode(X.hist.combine) <- mode(Y.hist.combine) <- mode(X.cur.con) <- mode(y.control) <- mode(X.cur.trt) <- "numeric"
  

  # Fit models
  
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # MAP_ESS: target ESS for PS-MAP prior method
  # sigma: the fixed reference scale in ess() for approximating mixture of normal
  # prior.trt.var: hyperparameter for the normal prior on treatment mean
  
  if (outcome.type == 1) {
    regroup_data <- PS_regroup(
      X.hist.combine = X.hist.combine, Y.hist.combine = Y.hist.combine, 
      X.cur.con = X.cur.con, Y.cur.con = y.control, X.cur.trt = X.cur.trt, 
      p = p, S = S, type = type
    )
    
    PS_MAP_fit <- PS_MAP.fit2(
      ybar.trt = mean(current.trt$Y), 
      sigma.hat.trt = var(current.trt$Y) * (n.trt - 1) / n.trt, 
      n.trt = n.trt, 
      tau.init = 1, 
      target.ESS = MAP_ESS, 
      sigma = sd(hist.data$Y),
      n.cur = regroup_data$n.cur, ybar.hist = regroup_data$Ybar.hist,
      SE.hist = regroup_data$SE.hist, overlap_coefficients = regroup_data$rS,
      lim = limits, data.indiv = y.control,
      prior.trt.var = 1e6, 
      alpha.sigma = 1, beta.sigma = 1, niter = niter
    )
  } else if (outcome.type == 2) {
    regroup_data <- PS_regroup_bin(
      X.hist.combine = X.hist.combine, Y.hist.combine = Y.hist.combine,
      X.cur.con = X.cur.con, Y.cur.con = y.control, X.cur.trt = X.cur.trt,
      p = p, S = S, type = type
    )
    
    PS_MAP_fit <- PS_MAP.fit2.bin(
      ybar.trt = sum(current.trt$Y),
      n.trt = n.trt,
      tau.init = 1,
      target.ESS = MAP_ESS,
      n.cur = regroup_data$n.cur,
      ybar.hist = regroup_data$Ybar.hist,
      n.hist = regroup_data$n.hist,
      overlap_coefficients = regroup_data$rS,
      lim = limits, data.indiv = y.control,
      alpha0 = 1, beta0 = 1, niter = niter
    )
  }
  
  
  PS_MAP_model <- PS_MAP_fit[[1]]
  PS_MAP_ESS <- PS_MAP_fit[[2]]
  
  # Gather Results
  
  results <- list(PSMAP = NULL, ESS = NULL)
  # results=list(PSMAP=NULL,MAP=NULL,ESS=NULL)
  
  results$PSMAP <- PS_MAP_model
  results$ESS <- PS_MAP_ESS
  
  if (getlimits == TRUE) {
    limits <- c(PS_MAP_fit[[3]] - 0.25, PS_MAP_fit[[3]] + 0.25)
    results <- limits
  }
  
  return(results)
  
}

ATE.psmap <- function(hist.data, current.ctrl, current.trt, MAP_ESS, S, 
                      type, outcome.type, showESS = FALSE){
  
  # get the approximated range for tau.scale in PS-MAP method
  limits <- run.psmap(
    hist.data = hist.data, current.ctrl = current.ctrl, current.trt = current.trt, 
    MAP_ESS = MAP_ESS, outcome.type = outcome.type, S = S, getlimits = TRUE, type = type
  )
  
  results <- run.psmap(
    hist.data = hist.data, current.ctrl = current.ctrl, current.trt = current.trt, 
    MAP_ESS = MAP_ESS, outcome.type = outcome.type, S = S, limits = limits, getlimits = FALSE, type = type
  )
  
  
  ATE_draws <- results[[1]][, "effect"]
  output <- c(
    ATE = mean(ATE_draws[[1]]), 
    prob = 2 * min(mean(ATE_draws[[1]] < 0), mean(ATE_draws[[1]] > 0))
    )
  
  if(showESS){
    output <- c(output, results$ESS)
  }
  
  return(output)
  
}


