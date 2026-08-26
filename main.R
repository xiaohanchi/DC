rm(list = ls())
# Mac: R packages (dbarts) and PyTorch each ship libomp; without this, import("ctgan") aborts the session
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE", OMP_NUM_THREADS = "1", PYTHONNOUSERSITE = "1")
library(pacman)
p_load(
  abind, arf, boot, bnlearn, broom, Brobdingnag, coda, cramer, dbarts, 
  distributional, dplyr, filelock, foreach, lmtest, MASS, MatchIt, 
  mvtnorm, sandwich, tableone, tibble, psych, purrr, tidyr, optmatch, 
  partitions, psrwe, randomForest, readr, ResourceSelection, rlist,
  rjags, runjags, RBesT, SAMprior, stringr, extraDistr, invgamma, modeest
)
select <- dplyr::select
rinvgamma <- invgamma::rinvgamma

### Settings =========
generative_model <- tibble(
  f.model = c(rep("bart", 3), rep("bnn", 3)),
  g.model = c("bn", "tvae", "arf", "bn", "tvae", "arf")
)
all.config <- expand.grid(
    model.type = c(2),
    bias.type = c(1),
    bias =  seq(-1, 1, 0.1),
    trt.eff = c(0, 0.5),
    scenarios =  22, #c(22, 24, 25, 27, 29, 30),
    sigma.rwdx = 1, 
    noise.rwd = 1,
    sigma.rctx = 1,
    sigma.rct = 1,
    rho.rwd = c(0.3), 
    exp.n = c(200), # c(100, 200),
    rwd.n = c(1000), # c(500, 1000),
    syn.nsample = 100, # c(100, 200, 500),
    syn.nset = 100,
    wt.type = 3, 
    wt.rho.y = 0.5,
    wt.b.y = c(2),  
    w0.val = -1, 
    bn.type = c(3),
    tvae_loss_fac =  1.25, 
    outcome.type = c(1), 
    var0.ess = c(0.05),
    prior.shrinkage = "dt(0, 5^(-2), 1)T(0,)",
    generative_idx = c(2, 3, 5, 6), 
    seed.pre = c(2344, 4566)
  )

### Setup =========

library(reticulate)
use_condaenv("r-reticulate", required = TRUE)
ctgan <- import("ctgan")

source("../data_generating.R")
source("../jags_functions.R")
source("../competing_functions.R")
source("../DC_functions.R")


### RUN code =======
resultpath.ATE <- "./results/output_ATE_sc00.csv"
resultpath.prob <- "./results/output_prob_sc00.csv"

for(rr in 1:2){
  seed <- as.numeric(paste0(all.config$seed.pre[sc00], seed00)) + c(0, 233)[rr]
  output.tmp <- MAIN.func(
    rwd.n = all.config$rwd.n[sc00], 
    exp.n = all.config$exp.n[sc00], 
    EHR.n = 2000, 
    synctrl.n = all.config$syn.nsample[sc00],
    trt.eff = all.config$trt.eff[sc00], 
    bias.c = all.config$bias[sc00], 
    syn.nset = all.config$syn.nset[sc00],
    scenario = all.config$scenarios[sc00], 
    var0.ess = all.config$var0.ess[sc00],
    prior.shrinkage = as.character(all.config$prior.shrinkage)[sc00],
    wt.type = all.config$wt.type[sc00], 
    wt.rho.y = all.config$wt.rho.y[sc00],
    wt.b.y = all.config$wt.b.y[sc00], 
    w0.val = all.config$w0.val[sc00], 
    rho.rwd = all.config$rho.rwd[sc00], 
    sigma.rwdx = all.config$sigma.rwdx[sc00], 
    sigma.rwd = all.config$noise.rwd[sc00], 
    sigma.rctx = all.config$sigma.rctx[sc00], 
    sigma.rct = all.config$sigma.rct[sc00],
    f.model = as.character(generative_model$f.model[all.config$generative_idx[sc00]]), 
    g.model = as.character(generative_model$g.model[all.config$generative_idx[sc00]]), 
    model.type = all.config$model.type[sc00], 
    bias.type = all.config$bias.type[sc00], 
    bn.type = all.config$bn.type[sc00], 
    tvae_loss_fac = all.config$tvae_loss_fac[sc00], 
    outcome.type = all.config$outcome.type[sc00],
    seed = seed, rep = (2 * rep00 - c(1, 0)[rr])
  )
  if (rr == 1) {
    output <- output.tmp
  } else {
    output <- bind_rows(output, output.tmp)
  }
}



### write file in a parallel way
lockfile <- "./results/lockfile.lock"
lock <- lock(lockfile, timeout = Inf)
if (!file.exists(resultpath.ATE)) {
  write_excel_csv(output$ATE, file = resultpath.ATE, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$ATE, file = resultpath.ATE, append = T, col_names = FALSE)
}

if (!file.exists(resultpath.prob)) {
  write_excel_csv(output$Prob, file = resultpath.prob, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$Prob, file = resultpath.prob, append = T, col_names = FALSE)
}
unlock(lock)
