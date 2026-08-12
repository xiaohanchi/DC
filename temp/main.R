rm(list = ls())
library(pacman)
p_load(
  BFI, dplyr, ggplot2, knitr, kableExtra, gridExtra, extraDistr, 
  filelock, readr, rjags, runjags, stringr
)
select <- dplyr::select

source("../jags_function.R")
source("../r_functions.R")


 
### RUN code ==========================
source("../scenarios.R")
all_config <- rbind(
  # continuous y
  # 4 sites w/ 5 periods
  expand.grid(
    y_type = 1,
    n_sc = c(1:2),
    beta_trt = c(0, 0.6), 
    site_delta_sc = c(1:3),
    drift_sc = c(1:3),
    active_time_sc = c(1),
    adjusted = c(TRUE, FALSE)
  ),
  # 6 sites w/ 8 periods
  expand.grid(
    y_type = 1,
    n_sc = c(4:5), 
    beta_trt = c(0, 0.5), 
    site_delta_sc = c(4:6),
    drift_sc = c(4:6),
    active_time_sc = c(2),
    adjusted = c(TRUE, FALSE)
  )
)

scenario <- list(
  beta0 = -1,
  beta = c(0.5, 0.3, -0.2),
  eps = 1,
  beta_trt = all_config$beta_trt[sc00],
  active_time = active_time[[all_config$active_time_sc[sc00]]],
  n_ctrl_by_site = n_by_site[[all_config$n_sc[sc00]]]$ctrl,
  n_trt_by_site = n_by_site[[all_config$n_sc[sc00]]]$trt,
  site_delta = site_delta[[all_config$site_delta_sc[sc00]]],
  drift = drift_mat[[all_config$drift_sc[sc00]]]
)

output <- main_func(
  pars = scenario, type = all_config$y_type[sc00], 
  adjusted = all_config$adjusted[sc00], lambda = 1/(10^2), 
  n_simu = 10, seed0 = seed00, rep = rep00
  )

all_res <- list(
  oneshotFP = output$beta_mat_FP,
  oneshotFP_noBorrow = output$beta_mat_FP_noBorrow,
  oneshotFP_indepTime = output$beta_mat_FP_indepTime,
  oneshotFP_noBorrow_indepTime = output$beta_mat_FP_noBorrow_indepTime,
  complete = output$beta_mat_complete,
  BFI = output$beta_mat_BFI,
  BFI_comp = output$beta_mat_BFI_comp,
  pooled = output$beta_mat_pool,
  local = output$beta_mat_local,
  localTM = output$beta_mat_localTM,
  poolTM = output$beta_mat_poolTM
)


### write file in a parallel way
lockfile <- "./results/lockfile.lock"
lock <- lock(lockfile, timeout = Inf)
for (method in names(all_res)) {
  resultpath <- file.path("./results", paste0("output_", method, "_sc", sc00, ".csv"))
  res_df <- as.data.frame(all_res[[method]])
  if (!file.exists(resultpath)) {
    write_excel_csv(res_df, file = resultpath, append = TRUE, col_names = TRUE)
  } else {
    write_excel_csv(res_df, file = resultpath, append = TRUE, col_names = FALSE)
  }
}
unlock(lock)