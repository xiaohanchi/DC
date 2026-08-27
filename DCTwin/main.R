rm(list = ls())
Sys.setenv(RETICULATE_PYTHON = "/rsrch4/home/biostatistics/xchi/.conda/envs/xchi/bin/python")
Sys.setenv(TABPFN_MODEL_CACHE_DIR = "/rsrch4/home/biostatistics/xchi/DCTwin/tabpfn-weights")
library(reticulate)
use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
reticulate::import("torch")
library(tabpfn)

library(pacman)
p_load(
  abind, boot, bnlearn, broom, Brobdingnag, coda, cramer, #dbarts, 
  distributional, dplyr, filelock, foreach, lmtest, MASS, MatchIt, 
  mvtnorm, sandwich, tableone, tibble, psych, purrr, tidyr, optmatch, 
  partitions, psrwe, randomForest, readr, ResourceSelection, rlist,
  rjags, runjags, RBesT, SAMprior, stringr, extraDistr, invgamma, modeest
)
select <- dplyr::select
rinvgamma <- invgamma::rinvgamma

source("../data_generating.R")
source("../jags_functions.R")
source("../competing_functions.R")
source("../DC_functions.R")

### Settings =========
all.config <- rbind(
  expand.grid(
    model.type = c(2, 4), # c(2, 4, 7),
    bias.type = c(1),
    bias =  0, # c(-0.6, -0.3, 0, 0.3, 0.6),
    trt.eff = c(0, 0.35),
    scenarios =  c(22, 24, 25, 27, 29, 30),
    sigma.rwdx = 1, 
    noise.rwd = 1,
    sigma.rctx = 1,
    sigma.rct = 1,
    rho.rwd = c(0.3), 
    exp.n = c(200), # c(100, 200),
    rwd.n = c(1000), # c(500, 1000),
    syn.nsample = 1000, # c(100, 200, 500),
    wt.type = 3, 
    wt.rho.x = c(10), 
    wt.b.x = c(1), 
    wt.rho.y = 0.5,
    wt.b.y = 2.5,  
    w0.val = -1,
    Gphi.type = c(3), #, c(1, 2, 3), 
    outcome.type = c(1), 
    var0.ess = c(0.05),
    seed.pre = c(1233, 2344, 3455, 4566, 5677)
  )
)
# all.config$trt.eff[all.config$trt.eff == 0.5] <- ifelse(all.config$exp.n[all.config$trt.eff == 0.5] == 100, 0.5, 0.35)
# Prior:  c("dunif(0, 10)", "dunif(0, 100)", "dt(0, 2.5^(-2), 1)T(0,)", "dt(0, 5^(-2), 1)T(0,)", "dt(0, 10^(-2), 1)T(0,)", "dt(0, 25^(-2), 1)T(0,)")

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
    scenario = all.config$scenarios[sc00], 
    var0.ess = all.config$var0.ess[sc00],
    wt.type = all.config$wt.type[sc00], 
    wt.rho.x = all.config$wt.rho.x[sc00],
    wt.b.x = all.config$wt.b.x[sc00], 
    wt.rho.y = all.config$wt.rho.y[sc00],
    wt.b.y = all.config$wt.b.y[sc00], 
    w0.val = all.config$w0.val[sc00], 
    rho.rwd = all.config$rho.rwd[sc00], 
    sigma.rwdx = all.config$sigma.rwdx[sc00], 
    sigma.rwd = all.config$noise.rwd[sc00], 
    sigma.rctx = all.config$sigma.rctx[sc00], 
    sigma.rct = all.config$sigma.rct[sc00],
    model.type = all.config$model.type[sc00], 
    bias.type = all.config$bias.type[sc00], 
    Gphi.type = all.config$Gphi.type[sc00], 
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
