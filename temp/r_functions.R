logit <- function(p) log(p / (1 - p))
expit <- function(x) exp(x) / (1 + exp(x))

make_xy <- function(dat, site_col = "site", y_col = "Y", intercept = TRUE) {
  y <- dat[[y_col]]
  x_cols <- setdiff(names(dat), c(site_col, y_col))
  X <- as.matrix(dat[, x_cols, drop = FALSE])
  if (intercept) X <- cbind(`(Intercept)` = 1, X)
  list(y = y, X = X, n = nrow(X), coef_names = colnames(X))
}


generate_data <- function(scenario = NULL, type = 1, 
                          adjusted = TRUE, seed = 233) {
  # type = 1 for continuous outcome; type = 2 for binary outcome
  # adjusted = TRUE for covariate adjusted model; adjusted = FALSE for unadjusted model
  if (is.null(scenario)) scenario <- get("scenario", envir = .GlobalEnv)

  pars <- scenario
  n_ctrl <- pars$n_ctrl_by_site
  n_trt <- pars$n_trt_by_site
  n_site <- length(n_ctrl)
  target_site <- which(n_trt > 0)

  if (length(n_trt) != n_site) stop("n_ctrl_by_site and n_trt_by_site must have the same length.")
  if (!is.null(seed)) set.seed(seed)

  data_list <- lapply(seq_len(n_site), function(kk) {
    n_k <- n_ctrl[kk] + n_trt[kk]
    trt_group <- c(rep(0, n_ctrl[kk]), rep(1, n_trt[kk])) %>% sample()
    ctrl_key <- paste0("site", kk, "_ctrl")
    trt_key <- paste0("site", kk, "_trt")
    temporal_ind <- if (ctrl_key %in% names(pars$active_time)) {
      tp <- rep(NA_integer_, n_k)
      tp[trt_group == 0] <- sample(rep(pars$active_time[[ctrl_key]], length.out = n_ctrl[kk]))
      tp[trt_group == 1] <- sample(rep(pars$active_time[[trt_key]], length.out = n_trt[kk]))
      tp
    } else {
      rep(pars$active_time[[paste0("site", kk)]], length.out = n_k) %>% sample()
    }
    X1 <- rnorm(n_k, mean = 0, sd = 2)
    X2 <- rnorm(n_k, mean = 2, sd = 3)
    X3 <- rbinom(n_k, size = 1, prob = 0.5)
    eta <- (pars$beta0 + pars$site_delta[kk]) + pars$drift[temporal_ind] +
      pars$beta[1] * X1 + pars$beta[2] * X2 + pars$beta[3] * X3 +
      pars$beta_trt * trt_group
    Y <- if (type == 1) {
      eta + rnorm(n_k, mean = 0, sd = pars$eps)
    } else {
      rbinom(n_k, size = 1, prob = expit(eta))
    }
    if(adjusted) {
      data <- tibble(site = kk, trt_group, temporal_ind, X1, X2, X3, Y) %>%
        arrange(trt_group, temporal_ind)
    } else {
      data <- tibble(site = kk, trt_group, temporal_ind, Y) %>%
        arrange(trt_group, temporal_ind)
    }
    
    data
  })

  list(
    data = bind_rows(data_list),
    scenario = pars,
    n_site = n_site,
    target_site = target_site,
    type = type
  )
}

run_MAP <- function(data, target_site, lambda = 0.0001, type) {
  # type: 1 = all sites, 2 = target site only
  family <- if (all(data$Y %in% c(0, 1))) "binomial" else "gaussian"

  if (type == 1) {
    M <- data %>% select(-c(site))
  } else if (type == 2) {
    M <- data %>%
      filter(site == target_site) %>%
      select(-c(site))
  }

  Lambda <- inv.prior.cov(
    as.data.frame(M %>% select(-c(Y))),
    lambda = lambda, family = family
  )
  fit <- MAP.estimation(
    M$Y,
    X = as.data.frame(M %>% select(-c(Y))), family = family, Lambda
  )

  return(fit)
}

run_TM <- function(data, target_site = NULL) {

  y_type <- if (all(data$Y %in% c(0, 1))) 2 else 1
  # target-site moments for continuous X
  ts_rows <- if (is.null(target_site)) rep(TRUE, nrow(data)) else data$site == target_site
  Xt <- data[ts_rows, setdiff(names(data), c("site", "Y", "temporal_ind", "trt_group")), drop = FALSE]
  Xt <- Xt[vapply(Xt, function(z) is.numeric(z) && length(unique(z)) > 2L, logical(1))]
  # standardize continuous X using target-site moments
  if (ncol(Xt)) {
    mu <- colSums(Xt) / nrow(Xt)
    s <- sqrt(pmax(colSums(Xt^2) / nrow(Xt) - mu^2, 1e-12))
    data[, names(mu)] <- scale(data[, names(mu), drop = FALSE], center = mu, scale = s)
  }
  X_mat <- data %>%
    select(-c(site, temporal_ind, Y)) %>%
    as.matrix()
  target_dat <- if (!is.null(target_site)) filter(data, site == target_site) else data
  X_target <- target_dat %>%
    filter(temporal_ind %in% intersect(temporal_ind[trt_group == 0], temporal_ind[trt_group == 1])) %>%
    select(-c(site, temporal_ind, Y)) %>%
    as.matrix()

  use_agg <- (y_type == 2 && ncol(X_mat) == 1L)
  if (y_type == 1) {
    # s_y <- sd(data$Y[if (is.null(target_site)) TRUE else data$site == target_site])
    jagsdata <- list(
      N = nrow(data),
      y = data$Y,
      Nperiod = max(data$temporal_ind),
      X = X_mat,
      time = data$temporal_ind,
      beta_p = ncol(X_mat),
      lambda_beta = 1 / 10^2
      # lambda_beta = (1 / 10^2) / max(s_y, 1e-8)^2
    )
  } else if (use_agg) {
    agg <- data %>%
      group_by(temporal_ind, trt_group) %>%
      summarise(y = sum(Y), n = dplyr::n(), .groups = "drop")
    jagsdata <- list(
      N_strata = nrow(agg),
      y = agg$y,
      n = agg$n,
      Nperiod = max(data$temporal_ind),
      time = agg$temporal_ind,
      trt = agg$trt_group
    )
  } else {
    jagsdata <- list(
      N = nrow(data),
      y = data$Y,
      Nperiod = max(data$temporal_ind),
      X = X_mat,
      time = data$temporal_ind,
      beta_p = ncol(X_mat),
      N_target = nrow(X_target),
      X_target = X_target
    )
  }
  data_glm <- data %>% select(-c(site))
  fit_glm <- tryCatch(
    glm(Y ~ .,
        data = transform(data_glm, temporal_ind = factor(temporal_ind)),
        family = if (y_type == 1) gaussian() else binomial()
    ),
    error = function(e) NULL
  )
  inits <- lapply(c(223, 323, 423, 523), function(s) {
    if (is.null(fit_glm)) {
      return(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    }
    cf <- coef(fit_glm)
    cf[!is.finite(cf)] <- 0
    out <- list(
      beta0 = unname(cf["(Intercept)"]),
      theta = replace(unname(cf[colnames(X_mat)]), is.na(cf[colnames(X_mat)]), 0),
      sd_alpha = 0.5,
      .RNG.name = "base::Mersenne-Twister", .RNG.seed = s
    )
    if (y_type == 1) out$log_sigma2 <- log(max(fit_glm$sigma^2, 1e-6))
    out
  })

  jagsmodel <- run.jags(
    model = if (y_type == 1) TM_continuous else if (use_agg) TM_binary_agg else TM_binary,
    monitor = if (y_type == 1) c("beta0", "theta", "prec_y") else c("beta0", "theta", "phat_ctrl", "phat_trt"),
    data = jagsdata, n.chains = 4,
    adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
    method = "rjags", plots = FALSE, silent.jags = T,
    inits = inits
  )

  if (y_type == 1) {
    theta_samples <- as.matrix(as.mcmc.list(jagsmodel, c("beta0", "theta", "prec_y")))
    colnames(theta_samples) <- c("beta0", colnames(X_mat), "prec_y")
    ATE.mcmc <- theta_samples[, "trt_group"]
  } else {
    theta_samples <- as.matrix(as.mcmc.list(jagsmodel, c("beta0", "theta", "phat_ctrl", "phat_trt")))
    colnames(theta_samples) <- c("beta0", colnames(X_mat), "p_ctrl", "p_trt")
    ATE.mcmc <- theta_samples[, "p_trt"] - theta_samples[, "p_ctrl"]
  }

  list(
    ATE = mean(ATE.mcmc),
    prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)),
    beta = colMeans(theta_samples[, 1:(1 + ncol(X_mat))]),
    cov = cov(theta_samples[, 1:(1 + ncol(X_mat)), drop = FALSE]),
    sigma2 = if (y_type == 1) 1 / (median(theta_samples[, "prec_y"])) else NULL
  )
}

run_localTM <- function(data) run_TM(data)
run_poolTM <- function(data, target_site) run_TM(data, target_site = target_site)


run_complete <- function(data, n_site, target_site) {

  y_type <- if (all(data$Y %in% c(0, 1))) 2 else 1
  # target-site moments for continuous X
  Xt <- data[data$site == target_site, setdiff(names(data), c("site", "Y", "temporal_ind", "trt_group")), drop = FALSE]
  Xt <- Xt[vapply(Xt, function(z) is.numeric(z) && length(unique(z)) > 2L, logical(1))]
  # standardize continuous X using target-site moments
  if (ncol(Xt)) {
    mu <- colSums(Xt) / nrow(Xt)
    s <- sqrt(pmax(colSums(Xt^2) / nrow(Xt) - mu^2, 1e-12))
    data[, names(mu)] <- scale(data[, names(mu), drop = FALSE], center = mu, scale = s)
  }
  X_mat <- data %>%
    select(-c(site, temporal_ind, Y)) %>%
    as.matrix()
  X_target <- data %>%
    filter(site == target_site) %>%
    filter(temporal_ind %in% intersect(temporal_ind[trt_group == 0], temporal_ind[trt_group == 1])) %>%
    select(-c(site, temporal_ind, Y)) %>%
    as.matrix()

  use_agg <- (y_type == 2 && ncol(X_mat) == 1L)
  if (y_type == 1) {
    # s_y <- sd(data$Y[data$site == target_site])
    jagsdata <- list(
      N = nrow(data),
      y = data$Y,
      Nperiod = max(data$temporal_ind),
      site = data$site,
      X = X_mat,
      time = data$temporal_ind,
      beta_p = ncol(X_mat),
      n_site = n_site,
      target_site = target_site,
      lambda_beta = 1 / 10^2
      # lambda_beta = (1 / 10^2) / max(s_y, 1e-8)^2
    )
  } else if (use_agg) {
    agg <- data %>%
      group_by(site, temporal_ind, trt_group) %>%
      summarise(y = sum(Y), n = dplyr::n(), .groups = "drop")
    jagsdata <- list(
      N_strata = nrow(agg),
      y = agg$y,
      n = agg$n,
      Nperiod = max(data$temporal_ind),
      site = agg$site,
      time = agg$temporal_ind,
      trt = agg$trt_group,
      n_site = n_site,
      target_site = target_site
    )
  } else {
    jagsdata <- list(
      N = nrow(data),
      y = data$Y,
      Nperiod = max(data$temporal_ind),
      site = data$site,
      X = X_mat,
      time = data$temporal_ind,
      beta_p = ncol(X_mat),
      n_site = n_site,
      target_site = target_site,
      N_target = nrow(X_target),
      X_target = X_target
    )
  }
  fit_glm <- tryCatch(
    glm(Y ~ .,
      data = transform(data, site = factor(site), temporal_ind = factor(temporal_ind)),
      family = if (y_type == 1) gaussian() else binomial()
    ),
    error = function(e) NULL
  )
  inits <- lapply(c(223, 323, 423, 523), function(s) {
    if (is.null(fit_glm)) {
      return(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
    }
    cf <- coef(fit_glm)
    cf[!is.finite(cf)] <- 0
    delta <- rep(0, n_site)
    for (nm in grep("^site\\d", names(cf), value = TRUE)) {
      delta[as.integer(sub("^site", "", nm))] <- cf[nm]
    }
    delta[target_site] <- 0
    out <- list(
      beta0 = unname(cf["(Intercept)"]),
      theta = replace(unname(cf[colnames(X_mat)]), is.na(cf[colnames(X_mat)]), 0),
      delta = delta, sd_alpha = 0.5, lambda = rep(1, n_site), tau = 1,
      .RNG.name = "base::Mersenne-Twister", .RNG.seed = s
    )
    if (y_type == 1) out$log_sigma2 <- log(max(fit_glm$sigma^2, 1e-6))
    out
  })

  jagsmodel <- run.jags(
    model = if (y_type == 1) complete_continuous else if (use_agg) complete_binary_agg else complete_binary,
    monitor = if (y_type == 1) c("beta0", "theta", "prec_y") else c("beta0", "theta", "phat_ctrl", "phat_trt"),
    data = jagsdata, n.chains = 4,
    adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
    method = "rjags", plots = FALSE, silent.jags = T,
    inits = inits
  )

  if (y_type == 1) {
    theta_samples <- as.matrix(as.mcmc.list(jagsmodel, c("beta0", "theta", "prec_y")))
    colnames(theta_samples) <- c("beta0", colnames(X_mat), "prec_y")
    ATE.mcmc <- theta_samples[, "trt_group"]
  } else {
    theta_samples <- as.matrix(as.mcmc.list(jagsmodel, c("beta0", "theta", "phat_ctrl", "phat_trt")))
    colnames(theta_samples) <- c("beta0", colnames(X_mat), "p_ctrl", "p_trt")
    ATE.mcmc <- theta_samples[, "p_trt"] - theta_samples[, "p_ctrl"]
  }

  result <- list(
    ATE = mean(ATE.mcmc),
    prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)),
    beta = colMeans(theta_samples[, 1:(1 + ncol(X_mat))]),
    cov = cov(theta_samples[, 1:(1 + ncol(X_mat)), drop = FALSE]),
    sigma2 = if (y_type == 1) 1 / (median(theta_samples[, "prec_y"])) else NULL
  )
  result
}

run_BFI <- function(data, n_site, target_site, lambda_loc = 0.0001, lambda_glb = 0.0001, homo = TRUE) {
  family <- if (all(data$Y %in% c(0, 1))) "binomial" else "gaussian"

  # target-site moments for continuous X
  Xt <- data[data$site == target_site, setdiff(names(data), c("site", "Y", "temporal_ind", "trt_group")), drop = FALSE]
  Xt <- Xt[vapply(Xt, function(z) is.numeric(z) && length(unique(z)) > 2L, logical(1))]
  target_x_summary <- list(n = nrow(Xt), sum = colSums(Xt), sumsq = colSums(Xt^2))

  Ms <- fits <- thetahats <- Ahats <- Lambdas <- list()
  warning_sites <- c()

  for (l in 1:n_site) {
    Ms[[l]] <- data %>%
      filter(site == l) %>%
      select(-c(site, temporal_ind))
    # standardize continuous X using target-site moments
    if (length(target_x_summary$sum)) {
      mu <- target_x_summary$sum / target_x_summary$n
      s <- sqrt(pmax(target_x_summary$sumsq / target_x_summary$n - mu^2, 1e-12))
      Ms[[l]][, names(mu)] <- scale(Ms[[l]][, names(mu), drop = FALSE], center = mu, scale = s)
    }
    Lambdas[[l]] <- inv.prior.cov(
      as.data.frame(Ms[[l]] %>% select(-c(Y))),
      lambda = lambda_loc, family = family
    )
    # unified prior: intercept/noise use same lambda as slopes
    # Lambdas[[l]]["(Intercept)", "(Intercept)"] <- lambda_loc / 100
    # if ("sigma2" %in% colnames(Lambdas[[l]])) Lambdas[[l]]["sigma2", "sigma2"] <- lambda_loc / 100
    fits[[l]] <- withCallingHandlers(
      MAP.estimation(
        y = Ms[[l]]$Y, X = as.data.frame(Ms[[l]] %>% select(-c(Y))),
        family = family, Lambda = Lambdas[[l]]
      ),
      warning = function(w) {
        warning_sites <<- c(warning_sites, l)
        invokeRestart("muffleWarning")
      }
    )
    thetahats[[l]] <- fits[[l]]$theta_hat
    Ahats[[l]] <- fits[[l]]$A_hat
  }

  if (homo) {
    Lambda_glb <- inv.prior.cov(
      X = as.data.frame(Ms[[1]] %>% select(-c(Y))),
      lambda = lambda_glb, family = family
    )
    # unified prior: intercept/noise use same lambda as slopes
    # Lambda_glb["(Intercept)", "(Intercept)"] <- lambda_glb / 100
    # if ("sigma2" %in% colnames(Lambda_glb)) Lambda_glb["sigma2", "sigma2"] <- lambda_glb / 100
    fit <- bfi(
      theta_hats = thetahats,
      A_hats = Ahats,
      Lambda = Lambda_glb,
      family = family,
      center_zero_sample = TRUE,
      which_cent_zeros = warning_sites,
      zero_sample_covs = "trt_group"
    )
  } else {
    Lambda_glb_hetero <- inv.prior.cov(
      X = as.data.frame(Ms[[1]] %>% select(-c(Y))),
      lambda = lambda_glb, family = family,
      stratified = TRUE, strat_par = 1, L = n_site
    )
    
    # i0 <- grep("^\\(Intercept\\)", colnames(Lambda_glb_hetero))
    # diag(Lambda_glb_hetero)[i0] <- lambda_glb / 100
    # if ("sigma2" %in% colnames(Lambda_glb_hetero)) Lambda_glb_hetero["sigma2", "sigma2"] <- lambda_glb / 100
    priors_all <- list(Lambdas[[1]], Lambda_glb_hetero)
    fit <- bfi(
      theta_hats = thetahats,
      A_hats = Ahats,
      Lambda = priors_all,
      family = family,
      center_zero_sample = TRUE,
      which_cent_zeros = warning_sites,
      zero_sample_covs = "trt_group",
      stratified = TRUE, strat_par = 1
    )
  }

  est <- if (homo) setNames(as.numeric(fit$theta_hat), colnames(fit$theta_hat)) else fit$theta_hat
  if (family == "gaussian") {
    fit$ATE <- unname(est["trt_group"])
    fit$prob <- 2 * (1 - pnorm(abs(fit$ATE / fit$sd["trt_group"])))
  } else if (family == "binomial") {
    beta_cols <- setdiff(names(Ms[[target_site]]), c("Y", "trt_group"))
    target_dat <- data %>%
      filter(site == target_site) %>%
      filter(temporal_ind %in% intersect(temporal_ind[trt_group == 0], temporal_ind[trt_group == 1]))
    X_target <- as.matrix(target_dat[, beta_cols, drop = FALSE])
    # standardize continuous X using target-site moments
    if (length(target_x_summary$sum)) {
      mu <- target_x_summary$sum / target_x_summary$n
      s <- sqrt(pmax(target_x_summary$sumsq / target_x_summary$n - mu^2, 1e-12))
      j <- intersect(names(mu), colnames(X_target))
      if (length(j)) X_target[, j] <- scale(X_target[, j, drop = FALSE], center = mu[j], scale = s[j])
    }
    intercept_col <- if (homo) "(Intercept)" else paste0("(Intercept)_loc", target_site)
    eta <- as.numeric(est[intercept_col] + X_target %*% est[beta_cols])
    p_est <- list(ctrl = expit(eta), trt = expit(eta + est["trt_group"]))
    fit$ATE <- mean(p_est$trt - p_est$ctrl)
    # delta method
    w <- p_est$trt * (1 - p_est$trt) - p_est$ctrl * (1 - p_est$ctrl)
    nm <- c(intercept_col, "trt_group", beta_cols)
    g <- c(
      mean(w), mean(p_est$trt * (1 - p_est$trt)), colMeans(w * X_target)
    )
    fit$prob <- 2 * (1 - pnorm(abs(fit$ATE / sqrt(c(t(g) %*% solve(fit$A_hat)[nm, nm] %*% g)))))
  }

  fit
}


# Surrogate log-likelihood for federated linear (continuous) or logistic (binary) regression.

run_surrogate <- function(
  data,
  data_site_i,
  beta_bar,
  sigma = NULL,
  site_col = "site",
  y_col = "Y",
  intercept = TRUE,
  maxit = 20
) {
  to_beta <- function(b, coef_names) {
    if (is.null(names(b))) as.numeric(b) else as.numeric(b[coef_names])
  }

  y_type <- if (all(data[[y_col]] %in% c(0, 1))) 2 else 1 # 1 = continuous, 2 = binary
  sites <- sort(unique(data[[site_col]]))
  local <- make_xy(data_site_i)
  bbar <- to_beta(beta_bar, local$coef_names)

  if (y_type == 1) {
    XtY_g <- XtX_g <- 0
    n_total <- rss <- 0

    for (s in sites) {
      prep <- make_xy(data[data[[site_col]] == s, ])
      XtY_g <- XtY_g + crossprod(prep$X, prep$y)
      XtX_g <- XtX_g + crossprod(prep$X)
      n_total <- n_total + prep$n
      b <- to_beta(beta_bar, prep$coef_names)
      rss <- rss + sum((prep$y - prep$X %*% b)^2)
    }

    if (is.null(sigma)) sigma <- sqrt(rss / nrow(data))
    inv_sigma2 <- 1 / sigma^2
    XtY_l <- crossprod(local$X, local$y)
    XtX_l <- crossprod(local$X)

    global_first <- as.vector(inv_sigma2 * (XtY_g - XtX_g %*% bbar) / n_total)
    global_second <- -inv_sigma2 * XtX_g / n_total
    local_first <- as.vector(inv_sigma2 * (XtY_l - XtX_l %*% bbar) / local$n)
    local_second <- -inv_sigma2 * XtX_l / local$n

    obj_fn <- function(beta) {
      b <- as.numeric(beta)
      local_ll <- sum((local$y - local$X %*% b)^2) / (2 * sigma^2 * local$n)
      linear <- sum((global_first - local_first) * b)
      quad <- as.vector(b - bbar) %*% (global_second - local_second) %*% as.vector(b - bbar)
      local_ll - linear - 0.5 * quad
    }
  } else {
    Xtgp_g <- XtWX_g <- 0
    n_total <- 0

    for (s in sites) {
      prep <- make_xy(data[data[[site_col]] == s, ])
      b <- to_beta(beta_bar, prep$coef_names)
      p <- as.vector(expit(prep$X %*% b))
      w <- p * (1 - p)
      Xtgp_g <- Xtgp_g + crossprod(prep$X, prep$y - p)
      XtWX_g <- XtWX_g + crossprod(prep$X, w * prep$X)
      n_total <- n_total + prep$n
    }

    p_l <- as.vector(expit(local$X %*% bbar))
    w_l <- p_l * (1 - p_l)

    global_first <- as.vector(Xtgp_g / n_total)
    global_second <- -XtWX_g / n_total
    local_first <- as.vector(crossprod(local$X, local$y - p_l) / local$n)
    local_second <- -crossprod(local$X, w_l * local$X) / local$n

    obj_fn <- function(beta) {
      b <- as.numeric(beta)
      p <- as.vector(expit(local$X %*% b))
      p <- pmax(pmin(p, 1 - 1e-15), 1e-15)
      local_ll <- -mean(local$y * log(p) + (1 - local$y) * log(1 - p))
      linear <- sum((global_first - local_first) * b)
      quad <- as.vector(b - bbar) %*% (global_second - local_second) %*% as.vector(b - bbar)
      local_ll - linear - 0.5 * quad
    }
  }

  est <- tryCatch(
    optim(bbar, obj_fn, method = "BFGS", control = list(maxit = maxit))$par,
    error = function(e) {
      message(e$message)
      rep(NA, length(bbar))
    }
  )

  setNames(est, local$coef_names)
}


run_oneshotFP <- function(data,
                          n_site,
                          site_col = "site",
                          y_col = "Y",
                          target_site = 4,
                          homo = TRUE,
                          homo_var = TRUE, 
                          no_borrow = TRUE,
                          rw_time = TRUE,
                          time_trend = TRUE,
                          lambda_local = 0.0001,
                          lambda_global = 0.0001,
                          sigma = NULL,
                          maxit = 50) {
  # homo: homo intercept
  # homo_var: homo model variance when y_type = 1 (continuous)
  
  Nperiod <- max(data$temporal_ind, na.rm = TRUE)
  y_type <- if (all(data[[y_col]] %in% c(0, 1))) 2 else 1 # 1 = continuous, 2 = binary
  y <- data[[y_col]]
  family <- if (all(y %in% c(0, 1))) "binomial" else "gaussian"

  if (time_trend) {
    x_cols <- setdiff(names(data), c(site_col, y_col))
    X_df <- cbind(
      as.data.frame(data[, setdiff(x_cols, "temporal_ind"), drop = FALSE]),
      setNames(
        as.data.frame(outer(data$temporal_ind, seq_len(Nperiod), `==`) + 0L),
        paste0("temporal_ind_", rev(seq_len(Nperiod)))
      )
    ) %>% select(-temporal_ind_1)
  } else {
    x_cols <- setdiff(names(data), c("temporal_ind", site_col, y_col))
    X_df <- as.data.frame(data[, x_cols, drop = FALSE])
  }
  coef_names <- c("(Intercept)", names(X_df))
  # shared <- if (y_type == 1) c(coef_names[-1], "sigma2") else coef_names[-1]
  shared <- coef_names[-1]

  Lambda_loc <- inv.prior.cov(
    X_df,
    lambda = lambda_local, family = family, intercept = TRUE
  )
  Lambda_global <- inv.prior.cov(
    X_df,
    lambda = lambda_global, family = family, intercept = TRUE
  )
  Lambda_glb_hetero <- inv.prior.cov(
    X_df,
    lambda = lambda_global, family = family,
    stratified = TRUE, strat_par = 1, L = n_site
  )
  
  Lambda_glb_hetero2 <- if (family == "gaussian") {
    inv.prior.cov(
      X_df,
      lambda = lambda_global, family = family,
      stratified = TRUE, strat_par = c(1, 2), L = n_site
    )
  } else {
    NULL
  }

  # target-site moments for continuous X
  Xt <- data[data[[site_col]] == target_site, setdiff(names(data), c(site_col, y_col, "temporal_ind", "trt_group")), drop = FALSE]
  Xt <- Xt[vapply(Xt, function(z) is.numeric(z) && length(unique(z)) > 2L, logical(1))]
  target_x_summary <- list(n = nrow(Xt), sum = colSums(Xt), sumsq = colSums(Xt^2))

  if (y_type == 1) {
    Hess_aa <- Hess_ab <- Hess_bb <- Hess_local <- eta_local <- vector("list", n_site)
    lambda_log_s2 <- lambda_local
    for (i in seq_len(n_site)) {
      if (time_trend) {
        prep <- make_xy(
          dat = data[data[[site_col]] == i, ],
          site_col = site_col, y_col = y_col, intercept = TRUE
        )
        X_df_i <- X_df[data[[site_col]] == i, , drop = FALSE]
        X_df_i <- cbind(`(Intercept)` = 1, X_df_i) %>% as.matrix()
      } else {
        prep <- make_xy(
          dat = data[data[[site_col]] == i, ] %>% select(-temporal_ind),
          site_col = site_col, y_col = y_col, intercept = TRUE
        )
        X_df_i <- prep$X
      }
      # standardize continuous X using target-site moments
      if (length(target_x_summary$sum)) {
        mu <- target_x_summary$sum / target_x_summary$n
        s <- sqrt(pmax(target_x_summary$sumsq / target_x_summary$n - mu^2, 1e-12))
        X_df_i[, names(mu)] <- scale(X_df_i[, names(mu), drop = FALSE], center = mu, scale = s)
      }

      p_coef <- ncol(X_df_i)
      beta_init <- lm.fit(X_df_i, prep$y)$coef
      beta_init[is.na(beta_init)] <- 0
      rss_init <- sum((prep$y - X_df_i %*% beta_init)^2)
      init_par <- c(beta_init, log_s2 = log(max(rss_init / prep$n, 1e-8)))
      fit <- optim(
        par = init_par,
        function(theta) {
          beta <- theta[seq_len(p_coef)]
          log_s2 <- theta[p_coef + 1]
          rss <- sum((prep$y - as.vector(X_df_i %*% beta))^2)
          nll <- rss / (2 * exp(log_s2)) + 0.5 * prep$n * log_s2
          prior_beta <- 0.5 * as.numeric(t(beta) %*% Lambda_loc[1:p_coef, 1:p_coef] %*% beta)
          prior_log_s2 <- 0.5 * lambda_log_s2 * log_s2^2
          nll + prior_beta + prior_log_s2
        },
        method = "L-BFGS",
        control = list(maxit = maxit)
      )
      beta_local_map <- fit$par[seq_len(p_coef)]
      log_s2_map <- fit$par[p_coef + 1]
      inv_s2_map <- exp(-log_s2_map)

      resids <- prep$y - X_df_i %*% beta_local_map
      hess_bb <- -inv_s2_map * crossprod(X_df_i) - Lambda_loc[1:p_coef, 1:p_coef]
      hess_bs <- -inv_s2_map * crossprod(X_df_i, resids)
      hess_ss <- -0.5 * inv_s2_map * sum(resids^2) - lambda_log_s2
      hess_logpost <- rbind(
        cbind(hess_bb, hess_bs), cbind(t(hess_bs), hess_ss)
      )

      Hess_bb[[i]] <- hess_logpost[1, 1, drop = FALSE]
      Hess_ab[[i]] <- hess_logpost[1, -1, drop = FALSE]
      Hess_aa[[i]] <- hess_logpost[-1, -1, drop = FALSE]
      Hess_local[[i]] <- hess_logpost
      eta_local[[i]] <- crossprod(hess_logpost, fit$par)
    }
    if (homo & homo_var) {
      Hess_global <- Reduce(`+`, Hess_local)
      eta_post <- Reduce(`+`, eta_local)
      hess_post <- Hess_global - Lambda_global + n_site * Lambda_loc
      beta_map <- solve(hess_post, eta_post)
      sigma2_map <- exp(beta_map[length(beta_map), 1])

      # current: homo is bundled with !time_trend model
      if (!time_trend) {
        result <- list(
          beta = beta_map[coef_names, 1],
          cov = solve(-hess_post),
          Hess = hess_post,
          Lambda = Lambda_global,
          Sigma_0 = solve(Lambda_global),
          sigma2 = sigma2_map
        )
      }
    } else if (!homo_var) {
      # hetero intercept + hetero sigma2
      p_all <- nrow(Hess_local[[1]])
      idx_b <- 1
      idx_s <- p_all
      idx_a <- 2:(p_all - 1)
      Ha <- Reduce(`+`, lapply(Hess_local, function(H) H[idx_a, idx_a, drop = FALSE]))
      Hess_global <- rbind(
        do.call(cbind, c(
          list(Ha),
          lapply(Hess_local, function(H) H[idx_a, idx_s, drop = FALSE]),
          lapply(Hess_local, function(H) H[idx_a, idx_b, drop = FALSE])
        )),
        do.call(rbind, lapply(seq_len(n_site), function(i) {
          H <- Hess_local[[i]]
          c(H[idx_s, idx_a], rep(0, i - 1), H[idx_s, idx_s], rep(0, n_site - i),
            rep(0, i - 1), H[idx_s, idx_b], rep(0, n_site - i))
        })),
        do.call(rbind, lapply(seq_len(n_site), function(i) {
          H <- Hess_local[[i]]
          c(H[idx_b, idx_a], rep(0, i - 1), H[idx_b, idx_s], rep(0, n_site - i),
            rep(0, i - 1), H[idx_b, idx_b], rep(0, n_site - i))
        }))
      )
      s2_nms <- paste0("sigma2_loc", seq_len(n_site))
      b0_nms <- paste0("(Intercept)_loc", seq_len(n_site))
      ord <- c(coef_names[-1], s2_nms, b0_nms)
      dimnames(Hess_global) <- list(ord, ord)
      eta_post <- c(
        Reduce(`+`, lapply(eta_local, function(e) e[idx_a])),
        sapply(eta_local, `[[`, idx_s),
        sapply(eta_local, `[[`, idx_b)
      ) %>% setNames(ord)
      Lam <- Lambda_glb_hetero2[ord, ord]
      Sigma_0 <- solve(Lam)
      hess_post <- Hess_global - Lam
      hess_post[shared, shared] <- hess_post[shared, shared] + n_site * Lambda_loc[shared, shared]
      site_pars <- c(s2_nms, b0_nms)
      hess_post[site_pars, site_pars] <- hess_post[site_pars, site_pars] +
        diag(c(rep(lambda_log_s2, n_site), rep(Lambda_loc[1, 1], n_site)))
      beta_map <- solve(hess_post, eta_post)

      jagsdata <- list(
        Nperiod = Nperiod,
        beta_p = length(shared[!grepl("^temporal_ind_", shared)]),
        shared_p = length(shared),
        n_sigma = n_site,
        n_site = n_site,
        target_site = target_site,
        n_p = length(beta_map),
        y_laplace = as.numeric(beta_map),
        invSigma = -(hess_post + Lam),
        lambda_beta = lambda_global
      )

      jags_model_cont <- if (no_borrow && rw_time) {
        FP_hetero_continuous
      } else if (no_borrow && !rw_time) {
        FP_hetero_continuous_indepTime
      } else if (!no_borrow && rw_time) {
        FP_continuous
      } else {
        FP_continuous_indepTime
      }
      jagsmodel <- run.jags(
        model = jags_model_cont,
        monitor = if (no_borrow) c("theta", "beta0_loc") else c("theta", "delta"),
        data = jagsdata, n.chains = 4,
        adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
        method = "rjags", plots = FALSE, silent.jags = T,
        inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
      )
      theta_samples <- as.matrix(as.mcmc.list(jagsmodel, "theta"))
      colnames(theta_samples) <- ord
      ord_beta <- setdiff(ord, s2_nms)
      ATE.mcmc <- theta_samples[, "trt_group"]

      result <- list(
        ATE = mean(ATE.mcmc),
        prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)),
        beta = setNames(colMeans(theta_samples[, ord_beta, drop = FALSE]), ord_beta),
        cov = cov(theta_samples[, ord_beta, drop = FALSE]),
        Hess = hess_post,
        Lambda_global = Lam,
        Sigma_0 = Sigma_0,
        sigma2 = exp(colMeans(theta_samples[, s2_nms, drop = FALSE]))
      )
    } else {
      Hess_global <- rbind(
        do.call(cbind, c(list(Reduce(`+`, Hess_aa)), lapply(Hess_ab, t))),
        do.call(rbind, lapply(seq_len(n_site), function(i) {
          c(Hess_ab[[i]], rep(0, i - 1), Hess_bb[[i]], rep(0, n_site - i))
        }))
      )
      ord <- c(coef_names[-1], "sigma2", paste0("(Intercept)_loc", seq_len(n_site)))
      dimnames(Hess_global) <- list(ord, ord)
      eta_post <- c(Reduce(`+`, lapply(eta_local, function(e) e[-1])), sapply(eta_local, `[[`, 1)) %>% setNames(ord)
      Sigma_0 <- solve(Lambda_glb_hetero[ord, ord])
      hess_post <- Hess_global - Lambda_glb_hetero[ord, ord]
      hess_post[shared, shared] <- hess_post[shared, shared] + n_site * Lambda_loc[shared, shared]
      hess_post["sigma2", "sigma2"] <- hess_post["sigma2", "sigma2"] + n_site * lambda_log_s2
      b0_nms <- paste0("(Intercept)_loc", seq_len(n_site))
      hess_post[b0_nms, b0_nms] <- hess_post[b0_nms, b0_nms] + diag(Lambda_loc[1, 1], n_site)
      beta_map <- solve(hess_post, eta_post)
      sigma2_map <- exp(beta_map["sigma2"])

      jagsdata <- list(
        Nperiod = Nperiod,
        beta_p = length(shared[!grepl("^temporal_ind_", shared)]),
        shared_p = length(shared),
        n_sigma = 1,
        n_site = n_site,
        target_site = target_site,
        n_p = length(beta_map),
        y_laplace = as.numeric(beta_map),
        invSigma = -(hess_post + Lambda_glb_hetero[ord, ord]),
        lambda_beta = lambda_global
      )

      # 2x2: (no_borrow x rw_time)
      # no_borrow=FALSE, rw_time=TRUE  -> horseshoe + RW
      # no_borrow=FALSE, rw_time=FALSE -> horseshoe + indep time
      # no_borrow=TRUE,  rw_time=TRUE  -> separate intercepts + RW
      # no_borrow=TRUE,  rw_time=FALSE -> separate intercepts + indep time
      jags_model_cont <- if (no_borrow && rw_time) {
        FP_hetero_continuous
      } else if (no_borrow && !rw_time) {
        FP_hetero_continuous_indepTime
      } else if (!no_borrow && rw_time) {
        FP_continuous
      } else {
        FP_continuous_indepTime
      }
      jagsmodel <- run.jags(
        model = jags_model_cont,
        monitor = if (no_borrow) c("theta", "beta0_loc") else c("theta", "delta"),
        data = jagsdata, n.chains = 4,
        adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
        method = "rjags", plots = FALSE, silent.jags = T,
        inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
      )
      theta_samples <- as.matrix(as.mcmc.list(jagsmodel, "theta"))
      colnames(theta_samples) <- ord
      ord_beta <- ord[ord != "sigma2"]
      ATE.mcmc <- theta_samples[, "trt_group"]

      result <- list(
        ATE = mean(ATE.mcmc),
        prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)),
        beta = setNames(colMeans(theta_samples[, ord_beta, drop = FALSE]), ord_beta),
        cov = cov(theta_samples[, ord_beta, drop = FALSE]),
        Hess = hess_post,
        Lambda_global = Lambda_glb_hetero[ord, ord],
        Sigma_0 = Sigma_0,
        sigma2 = exp(median(theta_samples[, "sigma2"]))
      )
    }
  } else if (y_type == 2) {
    Hess_aa <- Hess_ab <- Hess_bb <- Hess_local <- eta_local <- vector("list", n_site)
    for (i in seq_len(n_site)) {
      if (time_trend) {
        prep <- make_xy(
          dat = data[data[[site_col]] == i, ],
          site_col = site_col, y_col = y_col, intercept = TRUE
        )
        X_df_i <- X_df[data[[site_col]] == i, , drop = FALSE]
        X_df_i <- cbind(`(Intercept)` = 1, X_df_i) %>% as.matrix()
      } else {
        prep <- make_xy(
          dat = data[data[[site_col]] == i, ] %>% select(-temporal_ind),
          site_col = site_col, y_col = y_col, intercept = TRUE
        )
        X_df_i <- prep$X
      }
      # standardize continuous X using target-site moments
      if (length(target_x_summary$sum)) {
        mu <- target_x_summary$sum / target_x_summary$n
        s <- sqrt(pmax(target_x_summary$sumsq / target_x_summary$n - mu^2, 1e-12))
        X_df_i[, names(mu)] <- scale(X_df_i[, names(mu), drop = FALSE], center = mu, scale = s)
      }

      if (i == target_site) {
        site_dat <- data[data[[site_col]] == i, ]
        X_target <- X_df_i[site_dat$temporal_ind %in% intersect(site_dat$temporal_ind[site_dat$trt_group == 0], site_dat$temporal_ind[site_dat$trt_group == 1]), -1, drop = FALSE]
      }
      p_coef <- ncol(X_df_i)
      beta_init <- glm.fit(X_df_i, prep$y, family = binomial())$coef
      beta_init[!is.finite(beta_init)] <- 0
      if (max(abs(beta_init)) > 25) beta_init <- rep(0, p_coef)
      fit <- optim(
        par = beta_init,
        function(b) {
          p <- expit(as.vector(X_df_i %*% b))
          p <- pmax(pmin(p, 1 - 1e-15), 1e-15)
          nll <- -sum(prep$y * log(p) + (1 - prep$y) * log(1 - p))
          nll + 0.5 * as.numeric(t(b) %*% Lambda_loc %*% b) # MAP estimate
        },
        method = "BFGS",
        control = list(maxit = maxit)
      )

      beta_local_map <- fit$par
      p <- as.vector(expit(X_df_i %*% beta_local_map))
      score_logpost <- crossprod(X_df_i, (prep$y - p)) - crossprod(Lambda_loc, beta_local_map) # should be zero
      hess_logpost <- -crossprod(X_df_i, p * (1 - p) * X_df_i) - Lambda_loc # second derivative of log posterior
      Hess_bb[[i]] <- hess_logpost[1, 1, drop = FALSE]
      Hess_ab[[i]] <- hess_logpost[1, -1, drop = FALSE]
      Hess_aa[[i]] <- hess_logpost[-1, -1, drop = FALSE]

      Hess_local[[i]] <- hess_logpost
      eta_local[[i]] <- crossprod(hess_logpost, beta_local_map)
    }

    if (homo) {
      Hess_global <- Reduce(`+`, Hess_local)
      eta_post <- Reduce(`+`, eta_local)
      Sigma_0 <- solve(Lambda_global)
      hess_post <- Hess_global - Lambda_global + n_site * Lambda_loc

      beta_map <- solve(hess_post, eta_post)
      result <- list(
        beta = setNames(as.vector(beta_map), coef_names),
        cov = solve(-hess_post),
        Hess = hess_post,
        Lambda_global = Lambda_global,
        Sigma_0 = Sigma_0
      )
    } else {
      Hess_global <- rbind(
        do.call(cbind, c(list(Reduce(`+`, Hess_aa)), lapply(Hess_ab, t))),
        do.call(rbind, lapply(seq_len(n_site), function(i) {
          c(Hess_ab[[i]], rep(0, i - 1), Hess_bb[[i]], rep(0, n_site - i))
        }))
      )
      ord <- c(coef_names[-1], paste0("(Intercept)_loc", seq_len(n_site)))
      dimnames(Hess_global) <- list(ord, ord)
      eta_post <- c(Reduce(`+`, lapply(eta_local, function(e) e[-1])), sapply(eta_local, `[[`, 1)) %>% setNames(ord)
      Sigma_0 <- solve(Lambda_glb_hetero[ord, ord])
      hess_post <- Hess_global - Lambda_glb_hetero[ord, ord]
      hess_post[shared, shared] <- hess_post[shared, shared] + n_site * Lambda_loc[shared, shared]
      hess_post[setdiff(ord, shared), setdiff(ord, shared)] <- hess_post[setdiff(ord, shared), setdiff(ord, shared)] + diag(Lambda_loc[1, 1], nrow = n_site, ncol = n_site)

      beta_map <- solve(hess_post, eta_post)

      jagsdata <- list(
        Nperiod = Nperiod,
        beta_p = length(shared[!grepl("^temporal_ind_", shared)]),
        shared_p = length(shared),
        n_site = n_site,
        target_site = target_site,
        n_p = length(beta_map),
        y_laplace = as.numeric(beta_map),
        N_target = nrow(X_target),
        X_target = X_target,
        invSigma = -(hess_post + Lambda_glb_hetero[ord, ord])
      )

      # 2x2: (no_borrow x rw_time)
      # no_borrow=FALSE, rw_time=TRUE  -> horseshoe + RW
      # no_borrow=FALSE, rw_time=FALSE -> horseshoe + indep time
      # no_borrow=TRUE,  rw_time=TRUE  -> separate intercepts + RW
      # no_borrow=TRUE,  rw_time=FALSE -> separate intercepts + indep time
      jags_model_bin <- if (no_borrow && rw_time) {
        FP_hetero_binary
      } else if (no_borrow && !rw_time) {
        FP_hetero_binary_indepTime
      } else if (!no_borrow && rw_time) {
        FP_binary
      } else {
        FP_binary_indepTime
      }
      jagsmodel <- run.jags(
        model = jags_model_bin,
        monitor = if (no_borrow) c("theta", "beta0_loc", "phat_ctrl", "phat_trt") else c("theta", "delta", "phat_ctrl", "phat_trt"),
        data = jagsdata, n.chains = 4,
        adapt = 1000, burnin = 4000, sample = 5000, summarise = FALSE, thin = 2,
        method = "rjags", plots = FALSE, silent.jags = T,
        inits = lapply((c(1:4) * 100 + 123), function(s) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = s))
      )
      theta_samples <- as.matrix(as.mcmc.list(jagsmodel, c("theta", "phat_ctrl", "phat_trt")))
      colnames(theta_samples) <- c(ord, "p_ctrl", "p_trt")
      ord_beta <- ord
      ATE.mcmc <- theta_samples[, "p_trt"] - theta_samples[, "p_ctrl"]

      result <- list(
        ATE = mean(ATE.mcmc),
        prob = 2 * min(mean(ATE.mcmc < 0), mean(ATE.mcmc > 0)),
        beta = setNames(colMeans(theta_samples[, ord_beta, drop = FALSE]), ord_beta),
        cov = cov(theta_samples[, ord_beta, drop = FALSE]),
        Hess = hess_post,
        Lambda_global = Lambda_glb_hetero[ord, ord],
        Sigma_0 = Sigma_0
      )
    }
  }

  result
}


main_func <- function(
  pars = NULL,
  type = 1,
  adjusted = TRUE, 
  lambda = 0.0001,
  n_simu = 1000,
  seed0 = 233,
  rep = 1,
  verbose = TRUE
) {
  # type: = 1 for continuous outcome and = 2 for binary outcome

  runjags.options(silent.jags = TRUE, silent.runjags = TRUE, inits.warning = FALSE)

  target_site <- which(pars$n_trt_by_site > 0)
  beta0_true <- pars$beta0 + pars$site_delta[target_site]

  x_names <- if (is.null(names(pars$beta))) {
    paste0("X", seq_along(pars$beta))
  } else {
    names(pars$beta)
  }
  true_beta <- c(
    `(Intercept)` = beta0_true,
    trt_group = pars$beta_trt,
    # setNames(pars$beta, x_names)
    if (adjusted) setNames(pars$beta, x_names) else NULL
  )
  coef_names <- names(true_beta)
  shared_coef_names <- setdiff(coef_names, "(Intercept)")
  n_coef <- length(coef_names)
  col_names <- c(coef_names, "ATE", "prob")

  fill_hetero_beta <- function(mat, row, beta_vec) {
    mat[row, "(Intercept)"] <- beta_vec[paste0("(Intercept)_loc", target_site)]
    mat[row, shared_coef_names] <- beta_vec[shared_coef_names]
    mat
  }
  fill_pooled_beta <- function(mat, row, beta_vec) {
    mat[row, "(Intercept)"] <- beta_vec["beta0"]
    mat[row, shared_coef_names] <- beta_vec[shared_coef_names]
    mat
  }
  glm_to_hetero_beta <- function(fit, site_levels) {
    cf <- coef(fit)
    beta_vec <- cf[shared_coef_names]
    beta_vec[paste0("(Intercept)_loc", site_levels[1])] <- unname(cf["(Intercept)"])
    for (k in site_levels[-1]) {
      beta_vec[paste0("(Intercept)_loc", k)] <- unname(cf["(Intercept)"] + cf[paste0("site", k)])
    }
    beta_vec
  }

  beta_mat_FP <- beta_mat_FP_noBorrow <- beta_mat_FP_indepTime <- beta_mat_FP_noBorrow_indepTime <- beta_mat_FP_heteroVar <-
    beta_mat_BFI <- beta_mat_BFI_comp <- beta_mat_pool <- beta_mat_local <- beta_mat_localTM <- beta_mat_poolTM <- beta_mat_complete <-
    matrix(NA_real_, nrow = n_simu, ncol = length(col_names))
  colnames(beta_mat_FP) <- colnames(beta_mat_FP_noBorrow) <- colnames(beta_mat_FP_indepTime) <- colnames(beta_mat_FP_noBorrow_indepTime) <-
    colnames(beta_mat_FP_heteroVar) <-
    colnames(beta_mat_BFI) <- colnames(beta_mat_BFI_comp) <- colnames(beta_mat_pool) <- colnames(beta_mat_local) <- colnames(beta_mat_localTM) <- colnames(beta_mat_poolTM) <- colnames(beta_mat_complete) <- col_names

  for (r in seq_len(n_simu)) {
    if (verbose && r %% (n_simu / 10) == 0) message("Replicate ", r, " / ", n_simu)
    sim <- generate_data(
      scenario = pars, type = type, adjusted = adjusted, 
      seed = (seed0 * n_simu + r) * 10
      )

    lam_r <- lambda
    # lam_r <- if (type == 1) {
    #   s_y <- sd(sim$data$Y[sim$data$site == target_site[1]])
    #   lambda / max(s_y, 1e-8)^2
    # } else {
    #   lambda
    # }

    # proposed: Federated Platform Trial (borrow + RW time)
    fit_FP <- run_oneshotFP(
      data = sim$data, n_site = sim$n_site, target_site = target_site,
      homo = FALSE, homo_var = TRUE, no_borrow = FALSE, rw_time = TRUE, time_trend = TRUE,
      lambda_local = lam_r, lambda_global = lam_r
    )
    beta_mat_FP <- fill_hetero_beta(beta_mat_FP, r, fit_FP$beta)
    beta_mat_FP[r, c("ATE", "prob")] <- c(fit_FP$ATE, fit_FP$prob)

    # proposed with site-specific sigma2 (borrow + RW); binary has no sigma2 -> same as oneshotFP
    if (type == 1) {
      fit_FP_heteroVar <- run_oneshotFP(
        data = sim$data, n_site = sim$n_site, target_site = target_site,
        homo = FALSE, homo_var = FALSE, no_borrow = FALSE, rw_time = TRUE, time_trend = TRUE,
        lambda_local = lam_r, lambda_global = lam_r
      )
    } else {
      fit_FP_heteroVar <- fit_FP
    }
    beta_mat_FP_heteroVar <- fill_hetero_beta(beta_mat_FP_heteroVar, r, fit_FP_heteroVar$beta)
    beta_mat_FP_heteroVar[r, c("ATE", "prob")] <- c(fit_FP_heteroVar$ATE, fit_FP_heteroVar$prob)

    # proposed with separate intercepts (no borrowing)
    fit_FP_noBorrow <- run_oneshotFP(
      data = sim$data, n_site = sim$n_site, target_site = target_site,
      homo = FALSE, homo_var = TRUE, no_borrow = TRUE, rw_time = TRUE, time_trend = TRUE,
      lambda_local = lam_r, lambda_global = lam_r
    )
    beta_mat_FP_noBorrow <- fill_hetero_beta(beta_mat_FP_noBorrow, r, fit_FP_noBorrow$beta)
    beta_mat_FP_noBorrow[r, c("ATE", "prob")] <- c(fit_FP_noBorrow$ATE, fit_FP_noBorrow$prob)

    # borrow on intercept; independent time effects (no RW)
    fit_FP_indepTime <- run_oneshotFP(
      data = sim$data, n_site = sim$n_site, target_site = target_site,
      homo = FALSE, homo_var = TRUE, no_borrow = FALSE, rw_time = FALSE, time_trend = TRUE,
      lambda_local = lam_r, lambda_global = lam_r
    )
    beta_mat_FP_indepTime <- fill_hetero_beta(beta_mat_FP_indepTime, r, fit_FP_indepTime$beta)
    beta_mat_FP_indepTime[r, c("ATE", "prob")] <- c(fit_FP_indepTime$ATE, fit_FP_indepTime$prob)

    # no borrow on intercept or time (separate intercepts + indep time)
    fit_FP_noBorrow_indepTime <- run_oneshotFP(
      data = sim$data, n_site = sim$n_site, target_site = target_site,
      homo = FALSE, homo_var = TRUE, no_borrow = TRUE, rw_time = FALSE, time_trend = TRUE,
      lambda_local = lam_r, lambda_global = lam_r
    )
    beta_mat_FP_noBorrow_indepTime <- fill_hetero_beta(beta_mat_FP_noBorrow_indepTime, r, fit_FP_noBorrow_indepTime$beta)
    beta_mat_FP_noBorrow_indepTime[r, c("ATE", "prob")] <- c(fit_FP_noBorrow_indepTime$ATE, fit_FP_noBorrow_indepTime$prob)

    # complete data model
    fit_complete <- run_complete(
      data = sim$data, n_site = sim$n_site, target_site = target_site
    )
    beta_mat_complete <- fill_pooled_beta(beta_mat_complete, r, fit_complete$beta)
    beta_mat_complete[r, c("ATE", "prob")] <- c(fit_complete$ATE, fit_complete$prob)

    # BFI
    fit_BFI <- run_BFI(
      data = sim$data, n_site = sim$n_site, target_site = target_site,
      lambda_loc = lam_r, lambda_glb = lam_r,
      homo = FALSE
    )
    beta_mat_BFI <- fill_hetero_beta(beta_mat_BFI, r, fit_BFI$theta_hat)
    beta_mat_BFI[r, c("ATE", "prob")] <- c(fit_BFI$ATE, fit_BFI$prob)

    # target-site standardized data for MLE methods (Local / Pooled / BFI_comp)
    dat_std <- sim$data
    Xt <- dat_std[dat_std$site == target_site[1], setdiff(names(dat_std), c("site", "Y", "temporal_ind", "trt_group")), drop = FALSE]
    Xt <- Xt[vapply(Xt, function(z) is.numeric(z) && length(unique(z)) > 2L, logical(1))]
    if (ncol(Xt)) {
      mu <- colSums(Xt) / nrow(Xt)
      s <- sqrt(pmax(colSums(Xt^2) / nrow(Xt) - mu^2, 1e-12))
      dat_std[, names(mu)] <- scale(dat_std[, names(mu), drop = FALSE], center = mu, scale = s)
    }

    # BFI hetero intercept model with individual-level lm/glm (no time trend)
    data_BFI_comp <- dat_std %>%
      select(-temporal_ind) %>%
      mutate(site = factor(site))
    site_levels <- levels(data_BFI_comp$site)
    fit_formula <- reformulate(c(shared_coef_names, "site"), response = "Y")
    if (type == 1) {
      fit_BFI_comp <- lm(fit_formula, data = data_BFI_comp)
      beta_mat_BFI_comp <- fill_hetero_beta(beta_mat_BFI_comp, r, glm_to_hetero_beta(fit_BFI_comp, site_levels))
      beta_mat_BFI_comp[r, "ATE"] <- coef(fit_BFI_comp)["trt_group"]
      beta_mat_BFI_comp[r, "prob"] <- summary(fit_BFI_comp)$coefficients["trt_group", "Pr(>|t|)"]
    } else if (type == 2) {
      fit_BFI_comp <- glm(fit_formula, data = data_BFI_comp, family = "binomial")
      beta_mat_BFI_comp <- fill_hetero_beta(beta_mat_BFI_comp, r, glm_to_hetero_beta(fit_BFI_comp, site_levels))
      data_target <- dat_std %>%
        filter(site == target_site, temporal_ind %in% intersect(temporal_ind[trt_group == 0], temporal_ind[trt_group == 1])) %>%
        select(-temporal_ind) %>%
        mutate(site = factor(site, levels = site_levels))
      beta_mat_BFI_comp[r, "ATE"] <- mean(predict(fit_BFI_comp, newdata = transform(data_target, trt_group = 1), type = "response")) -
        mean(predict(fit_BFI_comp, newdata = transform(data_target, trt_group = 0), type = "response"))
      beta_mat_BFI_comp[r, "prob"] <- summary(fit_BFI_comp)$coefficients["trt_group", "Pr(>|z|)"]
    }

    # pooled
    data_pool <- dat_std %>% select(-c(site, temporal_ind))
    if (type == 1) {
      fit_pool <- lm(Y ~ ., data = data_pool)
      beta_mat_pool[r, coef_names] <- coef(fit_pool)
      beta_mat_pool[r, "ATE"] <- coef(fit_pool)["trt_group"]
      beta_mat_pool[r, "prob"] <- summary(fit_pool)$coefficients["trt_group", "Pr(>|t|)"]
    } else if (type == 2) {
      fit_pool <- glm(Y ~ ., data = data_pool, family = "binomial")
      beta_mat_pool[r, coef_names] <- coef(fit_pool)
      beta_mat_pool[r, "ATE"] <- mean(predict(fit_pool, newdata = transform(data_pool, trt_group = 1), type = "response")) -
        mean(predict(fit_pool, newdata = transform(data_pool, trt_group = 0), type = "response"))
      beta_mat_pool[r, "prob"] <- summary(fit_pool)$coefficients["trt_group", "Pr(>|z|)"]
    }

    # local: concurrent trt/ctrl periods only
    data_local <- dat_std %>%
      filter(site == target_site, temporal_ind %in% intersect(temporal_ind[trt_group == 0], temporal_ind[trt_group == 1])) %>%
      select(-c(site, temporal_ind))
    if (type == 1) {
      fit_local <- lm(Y ~ ., data = data_local)
      beta_mat_local[r, coef_names] <- coef(fit_local)
      beta_mat_local[r, "ATE"] <- coef(fit_local)["trt_group"]
      beta_mat_local[r, "prob"] <- summary(fit_local)$coefficients["trt_group", "Pr(>|t|)"]
    } else if (type == 2) {
      fit_local <- glm(Y ~ ., data = data_local, family = "binomial")
      beta_mat_local[r, coef_names] <- coef(fit_local)
      beta_mat_local[r, "ATE"] <- mean(predict(fit_local, newdata = transform(data_local, trt_group = 1), type = "response")) -
        mean(predict(fit_local, newdata = transform(data_local, trt_group = 0), type = "response"))
      beta_mat_local[r, "prob"] <- summary(fit_local)$coefficients["trt_group", "Pr(>|z|)"]
    }

    # local: time machine
    fit_localTM <- run_TM(sim$data %>% filter(site == target_site))
    beta_mat_localTM <- fill_pooled_beta(beta_mat_localTM, r, fit_localTM$beta)
    beta_mat_localTM[r, c("ATE", "prob")] <- c(fit_localTM$ATE, fit_localTM$prob)

    # pooled: time machine
    fit_poolTM <- run_TM(sim$data, target_site = target_site)
    beta_mat_poolTM <- fill_pooled_beta(beta_mat_poolTM, r, fit_poolTM$beta)
    beta_mat_poolTM[r, c("ATE", "prob")] <- c(fit_poolTM$ATE, fit_poolTM$prob)
  }

  add_result_cols <- function(mat, method) {
    if (!adjusted) {
      mat <- cbind(
        mat[, setdiff(colnames(mat), c("ATE", "prob")), drop = FALSE],
        matrix(NA_real_, nrow(mat), length(x_names), dimnames = list(NULL, x_names)),
        mat[, c("ATE", "prob"), drop = FALSE]
      )
    }
    cbind(replicate = (seq_len(nrow(mat)) + n_simu * (rep - 1)), method = method, mat)
  }

  list(
    beta_mat_FP = add_result_cols(beta_mat_FP, "oneshotFP"),
    beta_mat_FP_heteroVar = add_result_cols(beta_mat_FP_heteroVar, "oneshotFP_heteroVar"),
    beta_mat_FP_noBorrow = add_result_cols(beta_mat_FP_noBorrow, "oneshotFP_noBorrow"),
    beta_mat_FP_indepTime = add_result_cols(beta_mat_FP_indepTime, "oneshotFP_indepTime"),
    beta_mat_FP_noBorrow_indepTime = add_result_cols(beta_mat_FP_noBorrow_indepTime, "oneshotFP_noBorrow_indepTime"),
    beta_mat_complete = add_result_cols(beta_mat_complete, "Complete"),
    beta_mat_BFI = add_result_cols(beta_mat_BFI, "BFI"),
    beta_mat_BFI_comp = add_result_cols(beta_mat_BFI_comp, "BFI_comp"),
    beta_mat_pool = add_result_cols(beta_mat_pool, "Pooled"),
    beta_mat_local = add_result_cols(beta_mat_local, "Local"),
    beta_mat_localTM = add_result_cols(beta_mat_localTM, "LocalTM"),
    beta_mat_poolTM = add_result_cols(beta_mat_poolTM, "PoolTM")
  )
}
