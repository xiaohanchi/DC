library(tidyr)
library(readr)
library(doParallel)
library(foreach)
library(PMCMRplus)

solve.level <- function(rho, eff, tox) {
  c <- min(
    (rho * sqrt(eff * (1 - eff) * tox * (1 - tox)) + eff * tox),
    tox, eff
  )
  d <- min((eff - c), (1 - tox))
  a <- min((tox - c), (1 - eff))
  b <- 1 + round((c - eff - tox), 3)
  res <- round(c(a, b, c, d), 3)
  return(res)
}

expand.prob <- function(n1, n2, a0, b0) {
  a1 <- 0:n1
  a2 <- 0:n2
  a.value.comb <- expand_grid(a.value1 = a1, a.value2 = a2)
  b.value.comb <- expand_grid(b.value1 = rev(a1), b.value2 = rev(a2))
  combinations <- data.frame(
    a.value.grp1 = (a.value.comb$a.value1 + a0),
    b.value.grp1 = (b.value.comb$b.value1 + b0),
    a.value.grp2 = (a.value.comb$a.value2 + a0),
    b.value.grp2 = (b.value.comb$b.value2 + b0),
    probability = 0
  )
  combinations$probability <- sapply(1:dim(combinations)[1], function(r) {
    integrand <- function(x) {
      integrate.func.BCI(x,
        a = combinations$a.value.grp1[r], b = combinations$b.value.grp1[r],
        a.obd = combinations$a.value.grp2[r], b.obd = combinations$b.value.grp2[r]
      )
    }
    res <- integrate(integrand, lower = 0, upper = 1 - 1e-8)$value
    return(res)
  })
  return(combinations)
}


integrate.func.BCI <- function(x, a, b, a.obd, b.obd) {
  dbeta(x, a.obd, b.obd) * pbeta(x, a, b)
}

find.stage1.cutoff <- function(tox.upper, eff.lower, n.stage1, Cs2 = 0.8, step) {
  stopping.eff <- (which(
    pbeta(eff.lower, seq(0.1, (n.stage1 + 0.1)), seq((n.stage1 + 0.1), 0.1)) >=  Cs2
  ) %>% max()) - 1
  stage1.stop.eff <- stopping.eff + step
  stage1.cutoff.eff <- pbeta(eff.lower, (stage1.stop.eff + 0.1), (n.stage1 - stage1.stop.eff + 0.1))

  stage1.stop.tox <- (which(
    (1 - pbeta(tox.upper, seq(0.1, (n.stage1 + 0.1)), seq((n.stage1 + 0.1), 0.1))) >= stage1.cutoff.eff
  ) %>% min()) - 1

  stage1.cutoff.tox <- 1 - pbeta(tox.upper, (stage1.stop.tox + 0.1), (n.stage1 - stage1.stop.tox + 0.1))

  return(
    c(Cs1.tox = stage1.cutoff.tox, Cs1.eff = stage1.cutoff.eff)
  )
}

# n.stage1 <- 10
# n.stage2 <- 28
# N.dose <- 2
# N.arm <- 3
# tox.prob <- c(0.19, 0.02, 0.09, 0.63)
# eff.prob <- c(0.05, 0.38, 0.34, 0.63)
# rho <- 0.2
# n.simu <- 5000
# u.score <- c(0, 60, 40, 100)
# tox.upper <- 0.35
# eff.lower <- 0.25
# Cs <- 0.85


run.trial <- function(n.stage1, n.stage2,
                      N.dose, N.arm = 4,
                      tox.prob, eff.prob, rho = 0.2, u.score,
                      tox.upper = 0.35, eff.lower = 0.25, cutoff.step = 1,
                      Cs2 = 0.80, Ce.1, Ce.2, n.simu, upfront = TRUE) {
  set.seed(08291)

  a0.s1 <- rep(0.05, 4)
  a0.s2 <- b0.s2 <- 0.1
  multi.prob <- sapply(1:(N.dose + N.arm - 1), function(r) {
    solve.level(rho, eff.prob[r], tox.prob[r])
  })
  true.utility <- u.score %*% multi.prob[, N.arm:(dim(multi.prob)[2])]

  stage1.cutoff <- find.stage1.cutoff(
    tox.upper = tox.upper, eff.lower = eff.lower, n.stage1 = n.stage1, 
    Cs2 = Cs2, step = cutoff.step
  )
  Cs1.tox <- stage1.cutoff["Cs1.tox"]
  Cs1.eff <- stage1.cutoff["Cs1.eff"]

  ### stage I =====
  pi.hat.s1 <- array(0, dim = c((N.dose + N.arm - 1), nrow(multi.prob), n.simu))
  multi.Y.s1 <- sapply(1:(N.dose + N.arm - 1), function(r) {
    rmultinom(n.simu, n.stage1, multi.prob[, r])
  })
  Y.s1 <- array(t(multi.Y.s1), dim = c((N.dose + N.arm - 1), nrow(multi.prob), n.simu))
  if (!upfront) {
    Y.s1[1:(N.arm - 1), , ] <- 0
  }
  Y.tox.s1 <- Y.s1[, 1, ] + Y.s1[, 3, ]
  Y.eff.s1 <- Y.s1[, 3, ] + Y.s1[, 4, ]
  for (i in 1:nrow(multi.prob)) {
    pi.hat.s1[, i, ] <- (a0.s1[i] + Y.s1[, i, ]) / (sum(a0.s1) + n.stage1)
  }
  if (!upfront) {
    pi.hat.s1[1:(N.arm - 1), , ] <- 0
  }
  Pr.tox <- 1 - pbeta(
    tox.upper,
    (a0.s1[1] + a0.s1[3] + Y.tox.s1[N.arm:(dim(multi.prob)[2]), ]),
    (a0.s1[2] + a0.s1[4] + n.stage1 - Y.tox.s1[N.arm:(dim(multi.prob)[2]), ])
  )
  eff.a.value <- round(
    (a0.s1[3] + a0.s1[4] + Y.eff.s1[N.arm:(dim(multi.prob)[2]), ]),
    digits = 5
  )
  eff.b.value <- round(
    (a0.s1[1] + a0.s1[2] + n.stage1 - Y.eff.s1[N.arm:(dim(multi.prob)[2]), ]),
    digits = 5
  )
  Pr.eff <- pbeta(eff.lower, eff.a.value, eff.b.value)

  A.set <- which(Pr.tox < Cs1.tox & Pr.eff < Cs1.eff)
  for (i in 1:nrow(multi.prob)) {
    pi.hat.s1[N.arm:(dim(multi.prob)[2]), i, ][-A.set] <- -1
  }
  utility <- t(sapply(N.arm:(dim(multi.prob)[2]), function(r) {
    u.score %*% pi.hat.s1[r, , ]
  }))
  j.ast1 <- apply(utility, 2, which.max)
  j.ast1[which(apply(utility, 2, max) < 0)] <- -1 # no OBD; -A.set: exclude doses either terminated (proc.s1=0) or inadmissible

  ### stage II =====

  trial.s2.idx <- which(j.ast1 > 0)
  n.simu.s2 <- length(trial.s2.idx)
  A.indicator <- rep(1, n.simu)
  A.indicator[-trial.s2.idx] <- 0 # stopped trial

  Y.tox.s2 <- sapply(trial.s2.idx, function(r) {
    Y.tox.s1[c(1:(N.arm - 1), (N.arm - 1 + j.ast1[r])), r]
  })
  Y.eff.s2 <- sapply(trial.s2.idx, function(r) {
    Y.eff.s1[c(1:(N.arm - 1), (N.arm - 1 + j.ast1[r])), r]
  })

  Y.tox.s2 <- Y.tox.s2 + sapply(trial.s2.idx, function(r) {
    rbinom(rep(1, N.arm), n.stage2, prob = tox.prob[c(1:(N.arm - 1), (N.arm - 1 + j.ast1[r]))])
  })
  Y.eff.s2 <- Y.eff.s2 + sapply(trial.s2.idx, function(r) {
    rbinom(rep(1, N.arm), n.stage2, prob = eff.prob[c(1:(N.arm - 1), (N.arm - 1 + j.ast1[r]))])
  })

  #### BOP2-OC
  Pr.tox <- 1 - pbeta(
    tox.upper,
    (a0.s2 + Y.tox.s2[N.arm, ]),
    (b0.s2 + (n.stage1 + n.stage2) - Y.tox.s2[N.arm, ])
  )

  eff.a.value <- round((a0.s2 + Y.eff.s2), digits = 5)
  eff.b.value <- round((b0.s2 + (n.stage1 + n.stage2) - Y.eff.s2), digits = 5)

  Pr.eff <- pbeta(eff.lower, eff.a.value[N.arm, ], eff.b.value[N.arm, ])
  A.indicator[trial.s2.idx] <- (Pr.tox < Cs2 & Pr.eff < Cs2) * 1


  integrate.cases <- expand.prob(
    n1 = ifelse(upfront, (n.stage1 + n.stage2), n.stage2),
    n2 = (n.stage1 + n.stage2),
    a0 = a0.s2, b0 = b0.s2
  )
  BCI <- matrix(NA, nrow = (N.arm - 1), ncol = n.simu.s2)
  for (l in 1:(N.arm - 1)) {
    BCI[l, ] <- sapply(1:n.simu.s2, function(r) {
      idx <- c(which(integrate.cases$a.value.grp1 == eff.a.value[l, r] &
        integrate.cases$a.value.grp2 == eff.a.value[N.arm, r]), 1)
      res <- integrate.cases$probability[idx[1]]
      return(res)
    })
  }

  # 20260624 update
  # BCI.all <- cbind(BCI, matrix(-1,
  #   nrow = (N.arm - 1),
  #   ncol = (n.simu - n.simu.s2)
  # ))
  BCI.all <- matrix(-1, nrow = (N.arm - 1), ncol = n.simu)
  BCI.all[, trial.s2.idx] <- BCI

  j.ast <- j.ast1[j.ast1 > 0]
  total.n.s1 <- matrix(n.stage1, nrow = (N.dose + N.arm - 1), ncol = n.simu)
  n.cont <- sapply(1:n.simu.s2, function(r) {
    total.n.s1[(N.arm - 1 + j.ast[r]), trial.s2.idx[r]] <- (n.stage1 + n.stage2)
    return(total.n.s1[, trial.s2.idx[r]])
  })
  n.cont[1:(N.arm - 1), ] <- (n.stage1 + n.stage2)
  n.stop <- matrix(n.stage1, nrow = (N.dose + N.arm - 1), ncol = (n.simu - n.simu.s2))
  sample.size <- cbind(n.cont, n.stop)
  if (!upfront) {
    sample.size[(1:(N.arm - 1)), ] <- sample.size[(1:(N.arm - 1)), ] - n.stage1
  }
  avg.sample.size <- mean(colSums(sample.size))

  true.OBD <- which(
    true.utility == max(
      true.utility[which(tox.prob[-(1:(N.arm - 1))] <= tox.upper &
        eff.prob[-(1:(N.arm - 1))] >= eff.lower)]
    )
  )
  sel.opt.g <- sapply(1:n.simu, function(r) j.ast1[r] %in% true.OBD) * 1
  sel.opt.g <- which(sel.opt.g != 0) # correct groups

  CSP <- length(sel.opt.g) / n.simu
  power <- length(which(BCI.all[1, ] > Ce.1 & A.indicator == 1)) / n.simu # power
  opt.power <- length(which(BCI.all[1, sel.opt.g] > Ce.1 & A.indicator[sel.opt.g] == 1)) / length(sel.opt.g) # specific power
  GP <- CSP * opt.power

  if (N.arm == 4) {
    SR <- length(which(BCI.all[1, ] > Ce.1 &
      BCI.all[2, ] > Ce.2 &
      BCI.all[3, ] > Ce.2 &
      A.indicator == 1)) / n.simu
    opt.SR <- length(which(BCI.all[1, sel.opt.g] > Ce.1 &
      BCI.all[2, sel.opt.g] > Ce.2 &
      BCI.all[3, sel.opt.g] > Ce.2 &
      A.indicator[sel.opt.g] == 1)) / length(sel.opt.g)
  } else if (N.arm == 3) {
    SR <- length(which(BCI.all[1, ] > Ce.1 &
      BCI.all[2, ] > Ce.2 &
      A.indicator == 1)) / n.simu
    opt.SR <- length(which(BCI.all[1, sel.opt.g] > Ce.1 &
      BCI.all[2, sel.opt.g] > Ce.2 &
      A.indicator[sel.opt.g] == 1)) / length(sel.opt.g)
  } else if (N.arm == 2) {
    SR <- power
    opt.SR <- opt.power
  }
  GSR <- CSP * opt.SR

  summary.tab <- data.frame(
    Sample.Size = avg.sample.size,
    Power = power * 100, GP = GP * 100,
    SR = SR * 100, GSR = GSR * 100, CSP = CSP * 100
  )

  output <- list(
    true.utility = true.utility, true.OBD = true.OBD,
    summary.table = summary.tab
  )

  return(output)
}

dunnett.findcutoff <- function(n.stage1, n.stage2, p, N.arm, N.dose,
                               typeI.level = 0.10, NSR.level = 0.05, 
                               n.simu, upfront = TRUE, seed = 233) {
  n.arms <- N.arm + N.dose - 1
  n.total <- if (upfront) {
    n.stage1 * n.arms + n.stage2 * N.arm
  } else {
    n.stage1 * N.dose + n.stage2 * N.arm
  }
  n.base <- n.total %/% n.arms
  n.extra <- n.total %% n.arms
  n.per.arm <- rep(n.base, n.arms)
  if (n.extra > 0) {
    n.per.arm[seq_len(n.extra)] <- n.base + 1
  }
  pval.ctrl <- pval.single <- c()
  for (i in 1:n.simu) {
    set.seed(seed + 10 * i)
    data.tmp <- lapply(seq_len(n.arms), function(r) rbinom(n.per.arm[r], size = 1, prob = p[r]))
    data <- unlist(data.tmp)
    arm.label.tmp <- lapply(seq_len(n.arms), function(r) rep(r, n.per.arm[r]))
    arm.label <- unlist(arm.label.tmp)
    # comb vs ctrl
    # H0: comb <= SOC; a lower p-value means a better chance of comb > SOC
    dunnett.ctrl <- dunnettTest(
      x = data[which(arm.label %in% c(1, N.arm:(N.arm + N.dose - 1)))],
      g = arm.label[which(arm.label %in% c(1, N.arm:(N.arm + N.dose - 1)))],
      alternative = "greater"
    )
    pval.ctrl[i] <- min(dunnett.ctrl$p.value)

    # comb vs single
    if (N.arm == 4) {
      dunnett.single1 <- dunnettTest(
        x = data[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )

      dunnett.single2 <- dunnettTest(
        x = data[which(arm.label %in% c(3, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(3, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )
      pval.single[i] <- max(min(dunnett.single1$p.value), min(dunnett.single2$p.value))
    } else if (N.arm == 3) {
      dunnett.single <- dunnettTest(
        x = data[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )
      pval.single[i] <- min(dunnett.single$p.value)
    } else {
      pval.single[i] <- -1
    }
  }

  cutoff.ctrl <- quantile(pval.ctrl, typeI.level, type = 1)
  poc.idx <- which(pval.ctrl <= cutoff.ctrl)
  cutoff.single <- quantile(pval.single[poc.idx], NSR.level / typeI.level)

  return(
    c(ctrl = cutoff.ctrl, single = cutoff.single)
  )
}

run.dunnett <- function(n.stage1, n.stage2,
                        N.dose, N.arm = 4, eff.prob,
                        pval.cutoff, true.OBD, n.simu, upfront = TRUE, seed = 233) {
  n.arms <- N.arm + N.dose - 1
  n.total <- if (upfront) {
    n.stage1 * n.arms + n.stage2 * N.arm
  } else {
    n.stage1 * N.dose + n.stage2 * N.arm
  }
  n.base <- n.total %/% n.arms
  n.extra <- n.total %% n.arms
  n.per.arm <- rep(n.base, n.arms)
  if (n.extra > 0) {
    n.per.arm[seq_len(n.extra)] <- n.base + 1
  }

  pval.ctrl <- pval.single <- opt.dose <- c()
  for (i in 1:n.simu) {
    set.seed(seed + 10 * i)
    data.tmp <- lapply(seq_len(n.arms), function(r) rbinom(n.per.arm[r], size = 1, prob = eff.prob[r]))
    data <- unlist(data.tmp)
    arm.label.tmp <- lapply(seq_len(n.arms), function(r) rep(r, n.per.arm[r]))
    arm.label <- unlist(arm.label.tmp)
    # comb vs ctrl
    # H0: comb <= SOC; a lower p-value means a better chance of comb > SOC
    dunnett.ctrl <- dunnettTest(
      x = data[which(arm.label %in% c(1, N.arm:(N.arm + N.dose - 1)))],
      g = arm.label[which(arm.label %in% c(1, N.arm:(N.arm + N.dose - 1)))],
      alternative = "greater"
    )
    pval.ctrl[i] <- min(dunnett.ctrl$p.value)
    opt.dose[i] <- sapply(
      N.arm:(N.arm + N.dose - 1),
      function(r) sum(data[which(arm.label == r)])/n.per.arm[r]
    ) %>% which.max()

    # comb vs single
    if (N.arm == 4) {
      dunnett.single1 <- dunnettTest(
        x = data[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )

      dunnett.single2 <- dunnettTest(
        x = data[which(arm.label %in% c(3, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(3, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )
      pval.single[i] <- max(min(dunnett.single1$p.value), min(dunnett.single2$p.value))
    } else if (N.arm == 3) {
      dunnett.single <- dunnettTest(
        x = data[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        g = arm.label[which(arm.label %in% c(2, N.arm:(N.arm + N.dose - 1)))],
        alternative = "greater"
      )
      pval.single[i] <- min(dunnett.single$p.value)
    } else {
      pval.single[i] <- -1
    }
  }

  power <- length(which(pval.ctrl <= pval.cutoff[1])) / n.simu # power
  SR <- length(which(pval.ctrl <= pval.cutoff[1] & pval.single <= pval.cutoff[2])) / n.simu # power

  CSP <- length(which(opt.dose %in% true.OBD)) / n.simu
  GP <- length(which(pval.ctrl <= pval.cutoff[1] & opt.dose %in% true.OBD)) / n.simu
  GSR <- length(which(pval.ctrl <= pval.cutoff[1] &
    pval.single <= pval.cutoff[2] &
    opt.dose %in% true.OBD)) / n.simu # power

  summary.tab <- data.frame(
    Sample.Size = n.total,
    Power = power * 100, GP = GP * 100,
    SR = SR * 100, GSR = GSR * 100, CSP = CSP * 100
  )

  output <- list(
    summary.table = summary.tab
  )

  return(output)
}




# run.all <- function(Tox.prob, Eff.prob, n.scenarios) {
#   cl <- makeCluster(num.core - 1)
#   registerDoParallel(cl)
#   output <- foreach(
#     scenario.idx = 1:n.scenarios,
#     .export = c(
#       "run.trial", "solve.level", "expand.prob",
#       "integrate.func.BCI", "calibration.res"
#     ),
#     .packages = c("tidyr"), .combine = "rbind", .multicombine = TRUE
#   ) %dopar% {
#     res <- run.trial(
#       N.dose = 2, N.arm = 3,
#       n.stage1 = calibration.res$optimal$n.stage1,
#       n.stage2 = calibration.res$optimal$n.stage2,
#       tox.prob = as.numeric(Tox.prob[scenario.idx, ]),
#       eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
#       u.score = c(0, 60, 40, 100),
#       Ce.1 = calibration.res$optimal$Ce1,
#       Ce.2 = calibration.res$optimal$Ce2,
#       n.simu = 20000, upfront = TRUE
#     )
#     return(res$summary.table)
#   }
#   stopCluster(cl)
#   output <- round(output, digits = 2)
#   return(output)
# }

# RUN simulations ===========
## first get optimal cutoffs (exact values!)
# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/Projects/comb/code/calibration/sample_size/sample_size_cpp.R")
#
# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/Projects/comb/code/simulation/scenarios.R")
# num.core <- detectCores(logical = F)
#
# start <- Sys.time()
# results <- run.all(
#   Tox.prob = Tox.prob, Eff.prob = Eff.prob,
#   n.scenarios = n.scenarios
# )
# end <- Sys.time()
# end - start
