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
  n_simu = 100, seed0 = seed00, rep = rep00
  )

all_res <- list(
  oneshotFP = output$beta_mat_FP,
  complete = output$beta_mat_complete,
  BFI = output$beta_mat_BFI,
  combined = output$beta_mat_comb,
  local = output$beta_mat_local
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