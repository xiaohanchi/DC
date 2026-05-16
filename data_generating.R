logit <- function(p) log(p / (1 - p))
expit <- function(x) exp(x) / (1 + exp(x))

safe_factor <- function(data, col_names){
  cols <- intersect(col_names, names(data))
  data <- data %>% mutate(across(all_of(cols), factor))
  return(data)
}

bias.func <- function(X, c, type){
  switch (
    type,
    `1` = c, # constant
    `2` = X$X1 + c, # linear
    `3` = 0.5 * (X$X1 + X$X4) + c, 
    `4` = 0.5 * (X$X1 + X$X4) + X$X6 + c, 
    `5` = 0.3 * X$X1^2 + c, 
    `6` = 0.1 * X$X1^2 + 0.1 * exp(X$X4) + c,
    `7` = 0.1 * X$X1^2 + 0.1 * exp(X$X4) + X$X6 + c,
  )
}

generate.rwd <- function(N, trt.eff, bias.c, scenario, 
                         sigma.rwdx = 1, sigma.rwd = 1, rho = 0.3, 
                         model.type, bias.type, outcome.type = 1) {
  # sigma.rwdx: sd of X1-X3 in type 2 scenario 24
  # outcome.type: 1 for continuous; 2 for binary
  
  total.sc <- 30
  scenario <- (scenario - 1) %% total.sc + 1
  rho <- rho
  if (model.type == 1) {
    #### model type 1 ===========
    # linear, consistent relationship (7 covariates)
    if (scenario %in% c(1:12)) {
      covar.cont.lower <- covar.cont.upper <- matrix(NA, nrow = 4, ncol = 12)
      covar.cont.lower[, c(1:3, 6:8)] <- c(-2, -2, -2, -2)
      covar.cont.lower[, c(4, 5, 9, 10)] <- c(-2, -2, -4, -4)
      covar.cont.lower[, c(11, 12)] <- c(-4, -4, -2, -2)
      covar.cont.upper[, c(1:3, 4, 5, 9, 10)] <- c(2, 2, 2, 2)
      covar.cont.upper[, c(6:8)] <- c(2, 2, 4, 4)
      covar.cont.upper[, c(11, 12)] <- c(4, 4, 2, 2)

      covar.bin.prob <- matrix(NA, nrow = 3, ncol = 12)
      covar.bin.prob[, c(1:3, 4, 6, 8, 9, 11, 12)] <- c(0.5, 0.5, 0.5)
      covar.bin.prob[, c(5, 10)] <- c(0.2, 0.2, 0.5)
      covar.bin.prob[, c(7)] <- c(0.8, 0.8, 0.5)

      covar.cont <- sapply(1:4, function(r) {
        runif(n = N, min = covar.cont.lower[r, scenario], max = covar.cont.upper[r, scenario])
      })
      covar.bin <- sapply(1:3, function(r) rbinom(n = N, size = 1, prob = covar.bin.prob[r, scenario]))
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(
        bind_cols(covar.cont, covar.bin, covar.trt, covar.src,
          .name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)
        )
      )
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)

      coef.covar <- matrix(
        c(
          0.5, 0.5, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1
        ),
        ncol = 2, byrow = FALSE
      ) # quadratic, linear

      coef.trt <- trt.eff
      beta0 <- 2
      y.mean <- beta0 +
        (as.matrix(covar[, grep("X", colnames(covar))])^2) %*% matrix(coef.covar[, 1], ncol = 1) + # quad
        as.matrix(covar[, grep("X", colnames(covar))]) %*% matrix(coef.covar[, 2], ncol = 1) + # linear
        (covar.trt * coef.trt)
    } else if (scenario %in% c(13:21)) {
      if (scenario %in% c(13:19)) {
        marg.var <- rep(1, 7)
      } else if (scenario %in% c(20:21)) {
        marg.var <- c(2, 2, 1, 1, 1, 1, 1)^2
      }

      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var

      if (scenario %in% c(13:15, 20:21)) {
        covar.cont <- rmvnorm(n = N, mean = rep(0, 7), sigma = Covar.Mtx, method = "svd")
      } else if (scenario %in% c(16, 18)) {
        covar.cont <- rmvnorm(n = N, mean = c(0, 0, 1, 1, 0, 0, 0), sigma = Covar.Mtx, method = "svd")
      } else if (scenario %in% c(17, 19)) {
        covar.cont <- rmvnorm(n = N, mean = c(0, 0, 1, 1, 1, 1, 1), sigma = Covar.Mtx, method = "svd")
      }

      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)

      covar <- data.frame(cbind(covar.cont, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)),
        "treatment", "S"
      )
      covar <- tibble(covar)

      coef.covar <- matrix(
        c(
          0.5, 0.5, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1
        ),
        ncol = 2, byrow = FALSE
      )
      coef.trt <- trt.eff
      beta0 <- 2
      y.mean <- beta0 +
        (as.matrix(covar[, grep("X", colnames(covar))])^2) %*% matrix(coef.covar[, 1], ncol = 1) + # quad
        as.matrix(covar[, grep("X", colnames(covar))]) %*% matrix(coef.covar[, 2], ncol = 1) + # linear
        (covar.trt * coef.trt)
    } else if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5),
          rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2),
          rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)

      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
      coef.covar <- matrix(
        c(
          0.5, 0, 0, 0, 0, 0, 0,
          0.2, 0.5, 0.8, 0.2, 1, 1, 0.8
        ),
        ncol = 2, byrow = FALSE
      )
      coef.trt <- trt.eff
      beta0 <- 2
      y.mean <- beta0 +
        (as.matrix(covar[, grep("X", colnames(covar))])^2) %*% matrix(coef.covar[, 1], ncol = 1) + # quad
        as.matrix(covar[, grep("X", colnames(covar))]) %*% matrix(coef.covar[, 2], ncol = 1) + # linear
        (covar.trt * coef.trt)
    }
    
    if (outcome.type == 1) {
      y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
      Y.sd <- c()
      Y.sd[c(1, 4:7, 11, 22, 25, 27, 29:30)] <- 1
      Y.sd[c(2)] <- 3
      Y.sd[c(3, 8, 9, 10, 12, 26, 28)] <- 5
      Y.sd[c(13, 16, 17, 20)] <- 2
      Y.sd[c(14)] <- 4
      Y.sd[c(15, 18, 19, 21)] <- 6
      Y.sd[c(23, 24)] <- sigma.rwd
      data <- covar %>% mutate(Y = rnorm(
        n = N, mean = y.mean, sd = Y.sd[scenario]
      ))
    } else if (outcome.type == 2) {
      aa <- -1
      bb <- 0.2
      p.true <- expit(aa + bb * y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type))
      data <- covar %>% mutate(Y = rbern(n = N, prob = p.true))
    }
    
  } else if (model.type == 2) {
    ### model type 2 ===========
    # non-linear, consistent relationship (7 covariates)
    if (scenario %in% c(1:12)) {
      covar.cont.lower <- covar.cont.upper <- matrix(NA, nrow = 4, ncol = 12)
      covar.cont.lower[, c(1:3, 6:8)] <- c(-2, -2, -2, -2)
      covar.cont.lower[, c(4, 5, 9, 10)] <- c(-2, -2, -4, -4)
      covar.cont.lower[, c(11, 12)] <- c(-4, -4, -2, -2)
      covar.cont.upper[, c(1:3, 4, 5, 9, 10)] <- c(2, 2, 2, 2)
      covar.cont.upper[, c(6:8)] <- c(2, 2, 4, 4)
      covar.cont.upper[, c(11, 12)] <- c(4, 4, 2, 2)

      covar.bin.prob <- matrix(NA, nrow = 3, ncol = 12)
      covar.bin.prob[, c(1:3, 4, 6, 8, 9, 11, 12)] <- c(0.5, 0.5, 0.5)
      covar.bin.prob[, c(5, 10)] <- c(0.2, 0.2, 0.5)
      covar.bin.prob[, c(7)] <- c(0.8, 0.8, 0.5)

      covar.cont <- sapply(1:4, function(r) {
        runif(n = N, min = covar.cont.lower[r, scenario], max = covar.cont.upper[r, scenario])
      })
      covar.bin <- sapply(1:3, function(r) rbinom(n = N, size = 1, prob = covar.bin.prob[r, scenario]))
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)

      covar <- data.frame(
        bind_cols(covar.cont, covar.bin, covar.trt, covar.src,
          .name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)
        )
      )
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    } else if (scenario %in% c(13:21)) {
      if (scenario %in% c(13:19)) {
        marg.var <- c(1, 1, 1, 1, 0.5, 0.5, 0.5)^2
      } else if (scenario %in% c(20:21)) {
        marg.var <- c(2, 2, 1, 1, 0.5, 0.5, 0.5)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      if (scenario %in% c(13:15, 20:21)) {
        covar.cont <- rmvnorm(n = N, mean = rep(0, 7), sigma = Covar.Mtx, method = "svd")
      } else if (scenario %in% c(16, 18)) {
        covar.cont <- rmvnorm(n = N, mean = c(0, 0, 1, 1, 0, 0, 0), sigma = Covar.Mtx, method = "svd")
      } else if (scenario %in% c(17, 19)) {
        covar.cont <- rmvnorm(n = N, mean = c(0, 0, 1, 1, 1, 1, 1), sigma = Covar.Mtx, method = "svd")
      }

      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.cont, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    } else if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(sigma.rwdx, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }

    coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) -
      coef.covar[3] * covar$X3^2 - 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + # linear
      (covar.trt * coef.trt)
    
    
    if (outcome.type == 1) {
      y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
      Y.sd <- c()
      Y.sd[c(1, 4:7, 11, 13, 16, 17, 20, 22, 25, 27, 29:30)] <- 1 * sigma.rwd
      Y.sd[c(2, 14)] <- 1.5 * sigma.rwd
      Y.sd[c(3, 8, 9, 10, 12, 15, 18, 19, 21, 26, 28)] <- 2 * sigma.rwd
      Y.sd[c(23, 24)] <- 5 * sigma.rwd
      data <- covar %>% mutate(Y = rnorm(
        n = N, mean = y.mean, sd = Y.sd[scenario]
      ))
    } else if (outcome.type == 2) {
      aa <- -1.6
      bb <- 0.6
      p.true <- expit(aa + bb * y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type))
      data <- covar %>% mutate(Y = rbern(n = N, prob = p.true))
    }
    
    
  } else if (model.type == 3) {
    ### model type 3 ===========
    # non-linear, inconsistent relationship (7 covariates)
    if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) +
      coef.covar[3] * covar$X3^2 + 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + # linear
      (covar.trt * coef.trt)

    y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
    Y.sd <- rep(1, total.sc)
    Y.sd[c(22, 25, 27, 29:30)] <- 1
    Y.sd[c(26, 28)] <- 3
    Y.sd[c(23, 24)] <- sigma.rwd
    data <- covar %>% mutate(Y = rnorm(
      n = N, mean = y.mean, sd = Y.sd[scenario]
    ))
  } else if (model.type == 4) {
    ### model type 4 ===========
    # non-linear, consistent relationship (3/3 covariates)
    if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    
    # coef.covar <- c(2, 0, 0, 0, 0, 3)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 - 0.5 * covar$X3^2 + covar$X6 -  3 * covar$X6 * covar$X7 + (covar.trt * coef.trt)
    
    if (outcome.type == 1) {
      y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
      Y.sd <- c()
      Y.sd[c(1, 4:7, 11, 13, 16, 17, 20, 22, 25, 27, 29:30)] <- 1 * sigma.rwd
      Y.sd[c(2, 14)] <- 1.5 * sigma.rwd
      Y.sd[c(3, 8, 9, 10, 12, 15, 18, 19, 21, 26, 28)] <- 2 * sigma.rwd
      Y.sd[c(26, 28)] <- 5 * sigma.rwd
      Y.sd[c(23, 24)] <- 5 * sigma.rwd
      data <- covar %>% mutate(Y = rnorm(
        n = N, mean = y.mean, sd = Y.sd[scenario]
      )) %>% select(-c(X1, X2, X4, X5))
    } else if (outcome.type == 2) {
      aa <- -2
      bb <- 1
      p.true <- expit(aa + bb * y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type))
      data <- covar %>% mutate(Y = rbern(n = N, prob = p.true)) %>% 
        select(-c(X1, X2, X4, X5))
    }
    
  } else if (model.type == 5) {
    ### model type 5 ===========
    # non-linear, consistent relationship (3/7 covariates)
    if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    
    coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) +
      coef.covar[3] * covar$X3^2 + 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + 
      (covar.trt * coef.trt)
    
    y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
    Y.sd <- c()
    Y.sd[c(1, 4:7, 11, 13, 16, 17, 20, 22, 25, 27, 29:30)] <- 1
    Y.sd[c(2, 14)] <- 1.5
    Y.sd[c(3, 8, 9, 10, 12, 15, 18, 19, 21, 26, 28)] <- 2
    Y.sd[c(23, 24)] <- sigma.rwd
    data <- covar %>% mutate(Y = rnorm(
      n = N, mean = y.mean, sd = Y.sd[scenario]
    )) %>% select(-c(X2, X3, X4, X5))
  } else if (model.type == 6) {
    ### model type 6 ===========
    # non-linear, consistent relationship (7/3 covariates)
    if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    
    # coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 + 3 * expit(2 * covar$X1) + covar$X6 +  3 * covar$X6 * covar$X7 + (covar.trt * coef.trt)
    
    y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
    Y.sd <- c()
    Y.sd[c(1, 4:7, 11, 13, 16, 17, 20, 22, 25, 27, 29:30)] <- 1
    Y.sd[c(2, 14)] <- 1.5
    Y.sd[c(3, 8, 9, 10, 12, 15, 18, 19, 21, 26, 28)] <- 2
    Y.sd[c(23, 24)] <- sigma.rwd
    data <- covar %>% mutate(Y = rnorm(
      n = N, mean = y.mean, sd = Y.sd[scenario]
    ))
  } else if (model.type == 7) {
    ### model type 7 ===========
    # linear, consistent relationship (3/3 covariates)
    if (scenario %in% c(22:30)) {
      if (scenario %in% c(22:26, 29:30)) {
        marg.var <- rep(1, 3)^2
      } else if (scenario %in% c(27:28)) {
        marg.var <- c(2, 2, 2)^2
      }
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      # generate normal & uniform & binary covariates
      if (scenario %in% c(22:24, 27:30)) {
        covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
        )
      } else if (scenario %in% c(25, 26)) {
        covar.norm <- rmvnorm(n = N, mean = c(0, 0, 1), sigma = Covar.Mtx, method = "svd")
        covar.unif <- cbind(
          runif(n = N, min = -2, max = 4), runif(n = N, min = -1, max = 1)
        )
        covar.bin <- cbind(
          rbinom(n = N, size = 1, prob = 0.2), rbinom(n = N, size = 1, prob = 0.5)
        )
      }
      covar.trt <- rep(0, N)
      covar.src <- rep(1, N)
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
      coef.covar <- matrix(
        c(0, 0, 0, 0, 0, 0, 0,
          0, 0, 1.5, 0, 0, 1, 1), ncol = 2, byrow = FALSE
      )
      coef.trt <- trt.eff
      beta0 <- 2
      y.mean <- beta0 +
        (as.matrix(covar[, grep("X", colnames(covar))])^2) %*% matrix(coef.covar[, 1], ncol = 1) + # quad
        as.matrix(covar[, grep("X", colnames(covar))]) %*% matrix(coef.covar[, 2], ncol = 1) + # linear
        (covar.trt * coef.trt)
    }
    y.mean <- y.mean + bias.func(X = covar[grep("X", names(covar))], c = bias.c, type = bias.type)
    
    Y.sd <- c()
    Y.sd[c(22, 25, 27, 29:30)] <- 1 * sigma.rwd
    Y.sd[c(26, 28)] <- 5 * sigma.rwd
    Y.sd[c(23, 24)] <- 5 * sigma.rwd
    data <- covar %>% mutate(Y = rnorm(
      n = N, mean = y.mean, sd = Y.sd[scenario]
    )) %>% select(-c(X1, X2, X4, X5))
  }
  ## to keep the same seed for all scenarios
  data <- data %>%
    mutate(misX8 = rnorm(n = N, mean = 0, sd = 1), .after = X7) %>%
    mutate(misX9 = rnorm(n = N, mean = 0, sd = 1), .after = misX8)
  if (!(scenario %in% c(29))) {
    data <- data %>% select(-c(misX8, misX9))
  } 
  if (scenario %in% c(30)) {
    if(model.type %in% c(1:3, 6)) {
      data <- data %>% select(-c(X1, X4))
    } else if (model.type %in% c(4, 5, 7)) {
      data <- data %>% select(-c(X6, X7))
      }
  }
  return(data)
}


generate.RCT <- function(N, ratio, trt.eff, scenario, 
                         x.sd = 1, sd = 1, rho = 0.3, 
                         model.type, outcome.type = 1) {
  # ratio: randomization ratio; #trt/#total
  # ratio = 0 for all ctrl, and ratio = 1 for all trt
  # x.sd: sd in X1 - X3 in type 2
  # sd: noise in the data generating model
  # outcome.type: 1 for continuous; 2 for binary
  
  total.sc <- 30
  rho <- rho
  if (model.type == 1) {
    #### model type 1 ===========
    beta1q <- 0.5
    beta2q <- 0
    beta1 <- 0.2
    beta2 <- 0.5
    if ((scenario %% total.sc) %in% c(1:12)) {
      covar.cont <- sapply(1:4, function(r) runif(n = N, min = -2, max = 2))
      covar.bin <- sapply(1:3, function(r) rbinom(n = N, size = 1, prob = 0.5))
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.cont, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)

      coef.covar <- matrix(
        c(
          beta1q, beta2q, 0, 0, 0, 0, 0,
          beta1, beta1, 2, 1, 1, 3, 1
        ),
        ncol = 2, byrow = FALSE
      ) # quadratic, linear
      coef.trt <- trt.eff
    } else if ((scenario %% total.sc) %in% c(13:21)) {
      marg.var <- rep(1, 7)
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.cont <- rmvnorm(n = N, mean = rep(0, 7), sigma = Covar.Mtx)
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.cont, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
      coef.covar <- matrix(
        c(
          beta1q, beta2q, 0, 0, 0, 0, 0,
          beta1, beta1, 2, 1, 1, 3, 1
        ),
        ncol = 2, byrow = FALSE
      ) # quadratic, linear
      coef.trt <- trt.eff
    } else if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
      coef.covar <- matrix(
        c(
          beta1q, beta2q, 0, 0, 0, 0, 0,
          beta1, beta1, 0.8, 0.2, 1, 1, 0.8
        ),
        ncol = 2, byrow = FALSE
      ) # quadratic, linear
      coef.trt <- trt.eff
    }
    beta0 <- 2
    y.mean <- beta0 +
      (as.matrix(covar[, grep("X", colnames(covar))])^2) %*% matrix(coef.covar[, 1], ncol = 1) + # quad
      as.matrix(covar[, grep("X", colnames(covar))]) %*% matrix(coef.covar[, 2], ncol = 1) + # linear
      (covar.trt * coef.trt)

    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      data <- covar %>% mutate(Y = rnorm(
        n = N, mean = y.mean, sd = sd
      ))
    }
  } else if (model.type == 2) {
    ### model type 2 ===========
    if ((scenario %% total.sc) %in% c(1:12)) {
      covar.cont <- sapply(1:4, function(r) runif(n = N, min = -2, max = 2))
      covar.bin <- sapply(1:3, function(r) rbinom(n = N, size = 1, prob = 0.5))
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.cont, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    } else if ((scenario %% total.sc) %in% c(13:21)) {
      marg.var <- c(1, 1, 1, 1, 0.5, 0.5, 0.5)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.cont <- rmvnorm(n = N, mean = rep(0, 7), sigma = Covar.Mtx)
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.cont, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    } else if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(x.sd, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) -
      coef.covar[3] * covar$X3^2 - 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + # linear
      (covar.trt * coef.trt)
    
    if (outcome.type == 1) {
      data <- covar %>% mutate(Y = rnorm(n = N, mean = y.mean, sd = sd))
    } else if (outcome.type == 2) {
      aa <- -1.6
      bb <- 0.6
      p.true <- expit(aa + bb * y.mean)
      data <- covar %>% mutate(Y = rbern(n = N, prob = p.true))
    }
    
  } else if (model.type == 3) {
    ### model type 3 ===========
    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    coef.covar <- c(1, 1, 1.5, 0.1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) +
      coef.covar[3] * covar$X3^2 + 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + # linear
      (covar.trt * coef.trt)
    data <- covar %>% mutate(Y = rnorm(n = N, mean = y.mean, sd = sd))
  } else if (model.type == 4) {
    ### model type 4 ===========
    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    # coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 - 0.5 * covar$X3^2 + covar$X6 - 3 * covar$X6 * covar$X7 + (covar.trt * coef.trt)
    
    if (outcome.type == 1) {
      data <- covar %>% 
        mutate(Y = rnorm(n = N, mean = y.mean, sd = sd)) %>% 
        select(-c(X1, X2, X4, X5))
    } else if (outcome.type == 2) {
      aa <- -2
      bb <- 1
      p.true <- expit(aa + bb * y.mean)
      data <- covar %>% mutate(Y = rbern(n = N, prob = p.true)) %>% 
        select(-c(X1, X2, X4, X5))
    }
    
  } else if (model.type == 5) {
    ### model type 5 ===========
    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 +
      3 * expit(coef.covar[1] * covar$X1 + coef.covar[2] * covar$X2) +
      coef.covar[3] * covar$X3^2 + 0.3 * exp(coef.covar[4] * covar$X4) + coef.covar[5] * covar$X5 +
      coef.covar[6] * covar$X6 * covar$X7 + # linear
      (covar.trt * coef.trt)
    data <- covar %>% mutate(Y = rnorm(n = N, mean = y.mean, sd = sd)) %>% 
      select(-c(X2, X3, X4, X5))
  } else if (model.type == 6) {
    ### model type 6 ===========
    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    # coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 + 3 * expit(2 * covar$X1) + covar$X6 +  3 * covar$X6 * covar$X7 + (covar.trt * coef.trt)
    data <- covar %>% mutate(Y = rnorm(n = N, mean = y.mean, sd = sd))
  } else if (model.type == 7) {
    ### model type 7 ===========
    if ((scenario %% total.sc) %in% c(22:29, 0)) {
      marg.var <- rep(1, 3)^2
      Covar.Mtx <- rho * (sqrt(marg.var) %*% t(sqrt(marg.var)))
      diag(Covar.Mtx) <- marg.var
      covar.norm <- rmvnorm(n = N, mean = rep(0, 3), sigma = Covar.Mtx)
      covar.unif <- cbind(
        runif(n = N, min = -2, max = 2), runif(n = N, min = -1, max = 1)
      )
      covar.bin <- cbind(
        rbinom(n = N, size = 1, prob = 0.5), rbinom(n = N, size = 1, prob = 0.5)
      )
      covar.trt <- rep(0, N)
      covar.src <- rep(0, N)
      sample.idx <- sample(N, (N * ratio), replace = F)
      covar.trt[sample.idx] <- 1
      covar <- data.frame(cbind(covar.norm, covar.unif, covar.bin, covar.trt, covar.src))
      names(covar) <- c(
        sapply(1:(ncol(covar) - 2), function(r) paste0("X", r)), "treatment", "S"
      )
      covar <- tibble(covar)
    }
    # coef.covar <- c(1, 1, 0.5, 1, 1, -2)
    coef.trt <- trt.eff
    beta0 <- 2
    y.mean <- beta0 + 1.5 * covar$X3 + covar$X6 +  covar$X7 + (covar.trt * coef.trt)
    data <- covar %>% mutate(Y = rnorm(n = N, mean = y.mean, sd = sd)) %>% 
      select(-c(X1, X2, X4, X5))
  }

  data <- data %>%
    mutate(misX8 = rnorm(n = N, mean = 0, sd = 1), .after = X7) %>%
    mutate(misX9 = rnorm(n = N, mean = 0, sd = 1), .after = misX8)
  
  if (!((scenario %% total.sc) %in% c(29))) {
    data <- data %>% select(-c(misX8, misX9))
  } 
  if ((scenario %% total.sc) %in% c(0)) {
    if(model.type %in% c(1:3, 6)) {
      data <- data %>% select(-c(X1, X4))
    } else if (model.type %in% c(4, 5, 7)) {
      data <- data %>% select(-c(X6, X7))
    }
  }
  return(data)
}

prepare.data <- function(rwd.n, exp.n, EHR.n, trt.eff, bias.c, syn.nset, scenario, 
                         sigma.rwdx = 1, sigma.rwd = 1, sigma.rctx = 1, sigma.rct = 1, 
                         rho.rwd = 0.3, model.type, bias.type, outcome.type, seed) {
  total.sc <- 30
  set.seed(seed)
  rwd.data.raw <- generate.rwd(
    N = rwd.n, trt.eff = trt.eff, bias.c = bias.c, scenario = scenario, 
    sigma.rwdx = sigma.rwdx, sigma.rwd = sigma.rwd, rho = rho.rwd,
    model.type = model.type, bias.type = bias.type, outcome.type = outcome.type
  )
  # add EHR
  set.seed(seed + 2333)
  EHR.data <- generate.RCT(
    N = EHR.n, ratio = 0.5, trt.eff = trt.eff, scenario = scenario, 
    x.sd = sigma.rctx, sd = sigma.rct, rho = ifelse(rho.rwd == 0.3, rho.rwd, -rho.rwd), #20260214
    model.type = model.type, outcome.type = outcome.type
  ) 
  # RCT: a random subset from EHR?
  set.seed(seed + 6666)
  curr.data <- EHR.data %>% group_by(treatment) %>% 
    slice_sample(n = exp.n, replace = FALSE) %>% ungroup()
  true.ctrl.idx <- sample(which(curr.data$treatment == 0), size = (exp.n / 2), replace = FALSE)
  true.ctrl.s1 <- curr.data[true.ctrl.idx, ]
  exp.all <- curr.data %>% filter(treatment == 1)
  exp.all <- dplyr::select(exp.all, -c(S, treatment))
  true.ctrl.s1 <- dplyr::select(true.ctrl.s1, -c(S, treatment))
  curr.data <- curr.data[-true.ctrl.idx, ]
  
  ### true control from RCT 
  trueRCT <- list()
  for (ii in 1:syn.nset) {
    set.seed(seed + 10 * ii)
    ### generate true RCT for MAP
    trueRCT[[ii]] <- generate.RCT(
      N = (exp.n / 2), ratio = 0, trt.eff = trt.eff, scenario = scenario, 
      x.sd = sigma.rctx, sd = sigma.rct, 
      model.type = model.type, outcome.type = outcome.type
    ) %>% select(-c(treatment, S))
    
    if ((scenario %% total.sc) %in% c(1:12)) {
      trueRCT[[ii]] <- safe_factor(data = trueRCT[[ii]], col_names = c("X5", "X6", "X7"))
    } else if ((scenario %% total.sc) %in% c(22:29, 0)) {
      trueRCT[[ii]] <- safe_factor(data = trueRCT[[ii]], col_names = c("X6", "X7"))
    }
  }
  
  if ((scenario %% total.sc) %in% c(1:12)) {
    rwd.data.raw <- safe_factor(data = rwd.data.raw, col_names = c("X5", "X6", "X7"))
    exp.all <- safe_factor(data = exp.all, col_names = c("X5", "X6", "X7"))
    true.ctrl.s1 <- safe_factor(data = true.ctrl.s1, col_names = c("X5", "X6", "X7"))
    curr.data <- safe_factor(data = curr.data, col_names = c("X5", "X6", "X7"))
    EHR.data <- safe_factor(data = EHR.data, col_names = c("X5", "X6", "X7"))
  } else if ((scenario %% total.sc) %in% c(22:29, 0)) {
    rwd.data.raw <- safe_factor(data = rwd.data.raw, col_names = c("X6", "X7"))
    exp.all <- safe_factor(data = exp.all, col_names = c("X6", "X7"))
    true.ctrl.s1 <- safe_factor(data = true.ctrl.s1, col_names = c("X6", "X7"))
    curr.data <- safe_factor(data = curr.data, col_names = c("X6", "X7"))
    EHR.data <- safe_factor(data = EHR.data, col_names = c("X6", "X7"))
  }
  res <- list(
    rawRWD = rwd.data.raw,
    exp.all = exp.all,
    true.ctrl.s1 = true.ctrl.s1,
    EHR.data = EHR.data, 
    RCT.data = curr.data,
    trueRCT = trueRCT
  )
}

