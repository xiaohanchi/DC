rm(list = ls())
library(pacman)
p_load(
  abind, boot, bnlearn, broom, Brobdingnag, coda, cramer, dbarts, 
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
    methods = c("full", "selected")[2], # updated
    seed = seed
  )
  
  output <- list(
    ATE = mutate(res.s2$ATE, Replicate = rep, .before = 1),
    Prob = mutate(res.s2$Prob, Replicate = rep, .before = 1)
  )
  
  return(output)
}


### Settings =========
all.config <- rbind(
  expand.grid(
    model.type = c(2), # c(2, 4, 7),
    bias.type = c(1),
    bias =  seq(-0.6, 0.6, 0.1), # c(-0.6, -0.3, 0, 0.3, 0.6),
    trt.eff = c(0, 0.5),
    scenarios =  22,
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
    wt.rho.x = c(10), 
    wt.b.x = c(1), 
    wt.rho.y = 0.5,
    wt.b.y = c(2),  
    w0.val = -1, 
    bn.type = c(3), #, c(1, 2, 3), 
    outcome.type = c(1), 
    var0.ess = c(1e-10, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1),
    prior.shrinkage = "dt(0, 5^(-2), 1)T(0,)",
    seed.pre = c(2344, 4566)# c(1233, 2344, 3455, 4566, 5677)
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
    syn.nset = all.config$syn.nset[sc00],
    scenario = all.config$scenarios[sc00], 
    var0.ess = all.config$var0.ess[sc00],
    prior.shrinkage = as.character(all.config$prior.shrinkage)[sc00],
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
    bn.type = all.config$bn.type[sc00], 
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
