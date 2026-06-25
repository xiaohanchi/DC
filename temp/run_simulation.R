rm(list = ls())

all.config <- expand.grid(
  Cs2 = c(0.8, 0.85, 0.9, 0.95),
  cutoff.step = c(1, 2), 
  PowerMin = c(0.80, 0.85, 0.90), 
  SRMin = c(0.60, 0.65, 0.70)
)

Cs2 <- all.config$Cs2[idx00]
cutoff.step <- all.config$cutoff.step[idx00]
PowerMin <- all.config$PowerMin[idx00]
SRMin <- all.config$SRMin[idx00]

# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_research/Projects/comb/code/calibration/sample_size/sample_size_cpp.R")
source("~/comb2/calibration_sensitivity/sample_size_cpp.R")
# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_research/Projects/comb/code/simulation/run_trial.R")
source("~/comb2/calibration_sensitivity/run_trial.R")
num.core <- detectCores(logical = F)

#### upfront = TRUE ==========
##### case 1: N.arm = 4 ========

calibration.res <- all.search(
  n.stage1.lower = 6, n.stage1.upper = 40,
  n.stage2.lower = 2, n.stage2.upper = 40, step = 1,
  N.dose = 2, N.arm = 4, typeI.level = 0.10, NSR.level = 0.05,
  target.Power.min = PowerMin, target.SR.min = SRMin,
  eff.lower = 0.25, Cs2 = Cs2, cutoff.step = cutoff.step,
  diff.val = Inf, ratio.val = 0, upfront = TRUE
)
opt.para <- cbind(Case = 1, calibration.res$optimal)

pval.cutoff <- dunnett.findcutoff(
  n.stage1 = calibration.res$optimal$n.stage1,
  n.stage2 = calibration.res$optimal$n.stage2,
  p = rep(0.25, 5), N.arm = 4, N.dose = 2,
  typeI.level = 0.10, NSR.level = 0.05, 
  n.simu = 2e4, upfront = TRUE, seed = 233
)

# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_research/Projects/comb/code/simulation/scenarios_case1.R")
source("~/comb2/calibration_sensitivity/scenarios_case1.R")


run.all <- function(Tox.prob, Eff.prob, n.scenarios) {
  cl <- makeCluster(num.core - 1)
  registerDoParallel(cl)
  output <- foreach(
    scenario.idx = 1:n.scenarios,
    .export = c(
      "Cs2", "cutoff.step", "run.trial", "solve.level", "expand.prob",
      "integrate.func.BCI", "find.stage1.cutoff",
      "calibration.res", "pval.cutoff", "run.dunnett"
    ),
    .packages = c("tidyr", "PMCMRplus"), .combine = "rbind", .multicombine = TRUE
  ) %dopar% {
    res <- run.trial(
      N.dose = 2, N.arm = 4,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      tox.prob = as.numeric(Tox.prob[scenario.idx, ]),
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      u.score = c(0, 60, 40, 100), 
      tox.upper = 0.35, eff.lower = 0.25, 
      cutoff.step = cutoff.step, Cs2 = Cs2,
      Ce.1 = ceiling(calibration.res$optimal$Ce1 * 1e10) / 1e10,
      Ce.2 = ceiling(calibration.res$optimal$Ce2 * 1e10) / 1e10,
      n.simu = 20000, upfront = TRUE
    )

    res.dunnett <- run.dunnett(
      N.dose = 2, N.arm = 4,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      pval.cutoff = pval.cutoff, true.OBD = res$true.OBD, 
      n.simu = 20000, upfront = TRUE, seed = 233
    )

    return(
      rbind(res$summary.table, res.dunnett$summary.table)
      )
  }
  stopCluster(cl)
  output <- round(output, digits = 2)
  return(output)
}


results <- run.all(
  Tox.prob = Tox.prob, Eff.prob = Eff.prob,
  n.scenarios = n.scenarios
)

results.all <- cbind(Method = rep(c("BOP2-OC", "Dunnett"), n.scenarios), results)
# results.all %>% filter(Method == "BOP2-OC")

##### case 2 ========
calibration.res <- all.search(
  n.stage1.lower = 6, n.stage1.upper = 40,
  n.stage2.lower = 2, n.stage2.upper = 40, step = 1,
  N.dose = 2, N.arm = 3, typeI.level = 0.10, NSR.level = 0.05,
  target.Power.min = PowerMin, target.SR.min = SRMin,
  eff.lower = 0.25, Cs2 = Cs2, cutoff.step = cutoff.step,
  diff.val = Inf, ratio.val = 0, upfront = TRUE
)

opt.para <- rbind(opt.para, cbind(Case = 2, calibration.res$optimal))

pval.cutoff <- dunnett.findcutoff(
  n.stage1 = calibration.res$optimal$n.stage1,
  n.stage2 = calibration.res$optimal$n.stage2,
  p = rep(0.25, 4), N.arm = 3, N.dose = 2,
  typeI.level = 0.10, NSR.level = 0.05, n.simu = 2e4, upfront = TRUE
)

# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_research/Projects/comb/code/simulation/scenarios_case2.R")
source("~/comb2/calibration_sensitivity/scenarios_case2.R")

run.all <- function(Tox.prob, Eff.prob, n.scenarios) {
  cl <- makeCluster(num.core - 1)
  registerDoParallel(cl)
  output <- foreach(
    scenario.idx = 1:n.scenarios,
    .export = c(
      "Cs2", "cutoff.step", "run.trial", "solve.level", "expand.prob",
      "integrate.func.BCI", "find.stage1.cutoff",
      "calibration.res", "pval.cutoff", "run.dunnett"
    ),
    .packages = c("tidyr", "PMCMRplus"), .combine = "rbind", .multicombine = TRUE
  ) %dopar% {
    res <- run.trial(
      N.dose = 2, N.arm = 3,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      tox.prob = as.numeric(Tox.prob[scenario.idx, ]),
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      u.score = c(0, 60, 40, 100), tox.upper = 0.35,
      cutoff.step = cutoff.step, Cs2 = Cs2,
      Ce.1 = ceiling(calibration.res$optimal$Ce1 * 1e10) / 1e10,
      Ce.2 = ceiling(calibration.res$optimal$Ce2 * 1e10) / 1e10,
      n.simu = 20000, upfront = TRUE
    )
    
    res.dunnett <- run.dunnett(
      N.dose = 2, N.arm = 3,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      pval.cutoff = pval.cutoff, true.OBD = res$true.OBD, 
      n.simu = 20000, upfront = TRUE
    )
    
    return(
      rbind(res$summary.table, res.dunnett$summary.table)
    )
  }
  stopCluster(cl)
  output <- round(output, digits = 2)
  return(output)
}


results <- run.all(
  Tox.prob = Tox.prob, Eff.prob = Eff.prob,
  n.scenarios = n.scenarios
)

results.all <- cbind(results.all, NA, results)


#### upfront = FALSE ==========

##### case 1 ========

calibration.res <- all.search(
  n.stage1.lower = 6, n.stage1.upper = 40,
  n.stage2.lower = 2, n.stage2.upper = 40, step = 1,
  N.dose = 2, N.arm = 4, typeI.level = 0.10, NSR.level = 0.05,
  target.Power.min = PowerMin, target.SR.min = SRMin,
  eff.lower = 0.25, Cs2 = Cs2, cutoff.step = cutoff.step,
  diff.val = Inf, ratio.val = 0, upfront = FALSE
)

opt.para <- rbind(opt.para, cbind(Case = 3, calibration.res$optimal))

pval.cutoff <- dunnett.findcutoff(
  n.stage1 = calibration.res$optimal$n.stage1,
  n.stage2 = calibration.res$optimal$n.stage2,
  p = rep(0.25, 5), N.arm = 4, N.dose = 2,
  typeI.level = 0.10, NSR.level = 0.05, n.simu = 2e4, upfront = FALSE
)

# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_Research/Projects/comb/code/simulation/scenarios_case1.R")
source("~/comb2/calibration_sensitivity/scenarios_case1.R")

run.all <- function(Tox.prob, Eff.prob, n.scenarios) {
  cl <- makeCluster(num.core - 1)
  registerDoParallel(cl)
  output <- foreach(
    scenario.idx = 1:n.scenarios,
    .export = c(
      "Cs2", "cutoff.step", "run.trial", "solve.level", "expand.prob",
      "integrate.func.BCI", "find.stage1.cutoff",
      "calibration.res", "pval.cutoff", "run.dunnett"
    ),
    .packages = c("tidyr", "PMCMRplus"), .combine = "rbind", .multicombine = TRUE
  ) %dopar% {
    res <- run.trial(
      N.dose = 2, N.arm = 4,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      tox.prob = as.numeric(Tox.prob[scenario.idx, ]),
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      u.score = c(0, 60, 40, 100), tox.upper = 0.35,
      cutoff.step = cutoff.step, Cs2 = Cs2,
      Ce.1 = ceiling(calibration.res$optimal$Ce1 * 1e10) / 1e10,
      Ce.2 = ceiling(calibration.res$optimal$Ce2 * 1e10) / 1e10,
      n.simu = 20000, upfront = FALSE
    )
    
    res.dunnett <- run.dunnett(
      N.dose = 2, N.arm = 4,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      pval.cutoff = pval.cutoff, true.OBD = res$true.OBD, 
      n.simu = 20000, upfront = FALSE
    )
    
    return(
      rbind(res$summary.table, res.dunnett$summary.table)
    )
  }
  stopCluster(cl)
  output <- round(output, digits = 2)
  return(output)
}

results <- run.all(
  Tox.prob = Tox.prob, Eff.prob = Eff.prob,
  n.scenarios = n.scenarios
)


results.all <- cbind(results.all, NA, results)


##### case 2 ========

calibration.res <- all.search(
  n.stage1.lower = 6, n.stage1.upper = 40,
  n.stage2.lower = 2, n.stage2.upper = 40, step = 1,
  N.dose = 2, N.arm = 3, typeI.level = 0.10, NSR.level = 0.05,
  target.Power.min = PowerMin, target.SR.min = SRMin,
  eff.lower = 0.25, Cs2 = Cs2, cutoff.step = cutoff.step,
  diff.val = Inf, ratio.val = 0, upfront = FALSE
)

opt.para <- rbind(opt.para, cbind(Case = 4, calibration.res$optimal))

pval.cutoff <- dunnett.findcutoff(
  n.stage1 = calibration.res$optimal$n.stage1,
  n.stage2 = calibration.res$optimal$n.stage2,
  p = rep(0.25, 4), N.arm = 3, N.dose = 2,
  typeI.level = 0.10, NSR.level = 0.05, n.simu = 2e4, upfront = FALSE
)

# source("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_Research/Projects/comb/code/simulation/scenarios_case2.R")
source("~/comb2/calibration_sensitivity/scenarios_case2.R")

run.all <- function(Tox.prob, Eff.prob, n.scenarios) {
  cl <- makeCluster(num.core - 1)
  registerDoParallel(cl)
  output <- foreach(
    scenario.idx = 1:n.scenarios,
    .export = c(
      "Cs2", "cutoff.step", "run.trial", "solve.level", "expand.prob",
      "integrate.func.BCI", "find.stage1.cutoff",
      "calibration.res", "pval.cutoff", "run.dunnett"
    ),
    .packages = c("tidyr", "PMCMRplus"), .combine = "rbind", .multicombine = TRUE
  ) %dopar% {
    res <- run.trial(
      N.dose = 2, N.arm = 3,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      tox.prob = as.numeric(Tox.prob[scenario.idx, ]),
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      u.score = c(0, 60, 40, 100), tox.upper = 0.35,
      cutoff.step = cutoff.step, Cs2 = Cs2,
      Ce.1 = ceiling(calibration.res$optimal$Ce1 * 1e10) / 1e10,
      Ce.2 = ceiling(calibration.res$optimal$Ce2 * 1e10) / 1e10,
      n.simu = 20000, upfront = FALSE
    )
    
    res.dunnett <- run.dunnett(
      N.dose = 2, N.arm = 3,
      n.stage1 = calibration.res$optimal$n.stage1,
      n.stage2 = calibration.res$optimal$n.stage2,
      eff.prob = as.numeric(Eff.prob[scenario.idx, ]),
      pval.cutoff = pval.cutoff, true.OBD = res$true.OBD, 
      n.simu = 20000, upfront = FALSE
    )
    
    return(
      rbind(res$summary.table, res.dunnett$summary.table)
    )
  }
  stopCluster(cl)
  output <- round(output, digits = 2)
  return(output)
}

results <- run.all(
  Tox.prob = Tox.prob, Eff.prob = Eff.prob,
  n.scenarios = n.scenarios
)


results.all <- cbind(results.all, NA, results)


### Output Results ===========
opt.para <- cbind(
  Cs2 = Cs2, cutoff.step = cutoff.step, 
  Power.min = PowerMin, SR.min = SRMin, opt.para
  )
results.all <- cbind(
  Cs2 = Cs2, cutoff.step = cutoff.step, 
  Power.min = PowerMin, SR.min = SRMin, results.all
  )

parapath <- "./results/output_OptPara.csv"
resultspath <- "./results/output_OC.csv"

if (!file.exists(parapath)) {
  write_excel_csv(opt.para, file = parapath, na = " ", append = T, col_names = TRUE)
} else {
  write_excel_csv(opt.para, file = parapath, na = " ", append = T, col_names = FALSE)
}

if (!file.exists(resultspath)) {
  write_excel_csv(results.all, file = resultspath, na = " ", append = T, col_names = TRUE)
} else {
  write_excel_csv(results.all, file = resultspath, na = " ", append = T, col_names = FALSE)
}
