library(Rcpp)
library(dplyr)
library(tidyr)

localpath <- "~/Library/CloudStorage/OneDrive-InsideMDAnderson/PhD/01_research/Projects/comb/code/calibration/sample_size/DGP_comb_s1.cpp"
if(file.exists(localpath)){
  Rcpp::sourceCpp(localpath)
}else{
  Rcpp::sourceCpp("DGP_comb_s1.cpp")
}

one.search <- function(n.stage1, n.stage2,
                       N.dose, N.arm, eff.prob,
                       eff.lower, Cs2, cutoff.step,
                       typeI.level = 0.10, NSR.level = 0.05,
                       Ce.1 = 0, Ce.2 = 0, type = 1, find.cutoff = TRUE, upfront = TRUE) {
  # prob.inadm <- prob_inadm(
  #   Ndose = N.dose, Nstage1 = n.stage1, Nstage2 = n.stage2,
  #   prob_s1 = eff.prob[N.arm:length(eff.prob)], prob_s2 = max(eff.prob[N.arm:length(eff.prob)]),
  #   eff_lower = eff.lower, Cs2 = Cs2, cutoff_step = cutoff.step
  # )
  
  # type = 1 for max eff; = 2 for min eff
  temp <- generate_df(
    n_stage1 = n.stage1, n_stage2 = n.stage2,
    N_dose = N.dose, N_arm = N.arm,
    eff_prob = eff.prob, eff_lower = eff.lower,
    Cs2 = Cs2, cutoff_step = cutoff.step,
    type = type, upfront = upfront
  )

  all.cases.ctrl <- temp[, 1:8]
  quantile.df1 <- all.cases.ctrl[order(all.cases.ctrl$probability), ]
  quantile.df1$CDF <- sapply(1:dim(all.cases.ctrl)[1], function(r) sum(quantile.df1$PMF[1:r]))
  ##CDF: Prob(BCI<=x)

  if (isTRUE(find.cutoff)) {
    # 20260624
    # Ce.1 <- quantile.df1$probability[min(
    #   which(quantile.df1$CDF >= (max(quantile.df1$CDF) - typeI.level / (1 - prob.inadm)))
    # )]
    Ce.1 <- quantile.df1$probability[min(
      which(quantile.df1$CDF >= (max(quantile.df1$CDF) - typeI.level))
    )]
  }
  Power <- max(quantile.df1$CDF) - quantile.df1$CDF[min(which(quantile.df1$probability >= Ce.1))]
  quantile.df1 <- quantile.df1[quantile.df1$probability > Ce.1, ]
  sum.df1 <- quantile.df1 %>%
    group_by(a.value.grp2) %>%
    summarise(across(PMF, sum))
  colnames(sum.df1)[2] <- "PMF.grp2"


  if (N.arm != 2) {
    all.cases.single <- temp[, -(1:8)]
    colnames(all.cases.single) <- c(
      "a.value.grp1", "b.value.grp1", "a.value.grp2",
      "b.value.grp2", "probability", "PMF.grp1"
    )

    all.cases.single <- merge(all.cases.single, sum.df1, by = "a.value.grp2")
    all.cases.single$PMF <- all.cases.single$PMF.grp1 * all.cases.single$PMF.grp2

    quantile.df2 <- all.cases.single[order(all.cases.single$probability), ]
    quantile.df2$CDF <- sapply(1:dim(all.cases.single)[1], function(r) sum(quantile.df2$PMF[1:r]))


    if (isTRUE(find.cutoff)) {
      # 20260624
      # Ce.2 <- quantile.df2$probability[min(which(quantile.df2$CDF >= (Power - NSR.level / (1 - prob.inadm))))]
      Ce.2 <- quantile.df2$probability[min(which(quantile.df2$CDF >= (Power - NSR.level)))]
    }
    SR <- Power - quantile.df2$CDF[min(which(quantile.df2$probability >= Ce.2))]
  } else if (N.arm == 2) {
    if (isTRUE(find.cutoff)) {
      Ce.2 <- 0
    }
    SR <- Power
  }
  # 20260624
  # Power <- Power * (1 - prob.inadm)
  # SR <- SR * (1 - prob.inadm)

  res <- c(Ce.1 = Ce.1, Ce.2 = Ce.2, Power = Power, SR = SR)

  return(res)
}


# n.stage1.lower <- 2
# n.stage1.upper <- 30
# n.stage2.lower <- 2
# n.stage2.upper <- 30
# step <- 2
# N.dose <- 2
# N.arm <- 4
# typeI.level <- 0.10
# NSR.level <- 0.05
# diff.val <- Inf
# eff.prob <- eff.null

all.search <- function(n.stage1.lower, n.stage1.upper,
                       n.stage2.lower, n.stage2.upper, step = 2,
                       N.dose, N.arm = 4,
                       typeI.level = 0.10,
                       NSR.level = 0.05,
                       target.Power.min,
                       target.SR.min,
                       eff.lower, Cs2, cutoff.step = 1,
                       diff.val = Inf, ratio.val = 0, upfront = TRUE) {
  # diff.val: if requiring (n1 - n2) <= diff.val;
  # ratio.val: if requiring (n1 / n2) >= ratio.val;

  ###### setting

  eff.null <- rep(0.25, (N.dose + N.arm - 1))
  eff.alt <- c(0.25, rep(0.35, N.arm - 2), rep(0.55, N.dose))

  opt.idx <- TRUE

  sample.sizes <- expand.grid(
    n.stage1 = seq(n.stage1.lower, n.stage1.upper, step),
    n.stage2 = seq(n.stage2.lower, n.stage2.upper, step)
  )
  sample.sizes$diff <- sample.sizes$n.stage1 - sample.sizes$n.stage2
  sample.sizes <- dplyr::filter(sample.sizes, diff <= diff.val)
  sample.sizes <- sample.sizes[, -grep("diff", colnames(sample.sizes))]

  sample.sizes$ratio <- sample.sizes$n.stage1 / sample.sizes$n.stage2
  sample.sizes <- dplyr::filter(sample.sizes, ratio >= ratio.val)
  sample.sizes <- sample.sizes[, -grep("ratio", colnames(sample.sizes))]

  if (upfront) {
    sample.sizes$n.total.trial <- (N.dose + N.arm - 1) * sample.sizes$n.stage1 +
      N.arm * sample.sizes$n.stage2
  } else {
    sample.sizes$n.total.trial <- (N.dose) * sample.sizes$n.stage1 +
      N.arm * sample.sizes$n.stage2
  }
  sample.sizes <- sample.sizes[order(
    sample.sizes$n.total.trial,
    -sample.sizes$n.stage1
  ), ]
  rownames(sample.sizes) <- NULL
  if (dim(sample.sizes)[1] <= 200) {
    rough.idx <- quantile(1:dim(sample.sizes)[1], probs = seq(0, 1, 1 / 20), type = 1)
  } else {
    rough.idx <- seq(1, dim(sample.sizes)[1], 10)
  }

  # val <- 0
  # for (i in rough.idx) {
  #   val <- val + 1
  #   n.stage1 <- sample.sizes$n.stage1[i]
  #   n.stage2 <- sample.sizes$n.stage2[i]
  # 
  #   run.null.max <- one.search(
  #     n.stage1 = n.stage1, n.stage2 = n.stage2,
  #     N.dose = N.dose, N.arm = N.arm, eff.prob = eff.null,
  #     eff.lower = eff.lower, Cs2 = Cs2, cutoff.step = cutoff.step,
  #     typeI.level = typeI.level, NSR.level = NSR.level, type = 1, upfront = upfront
  #   )
  #   Ce.1 <- run.null.max["Ce.1"]
  #   Ce.2 <- run.null.max["Ce.2"]
  # 
  # 
  #   run.alt.min <- one.search(
  #     n.stage1 = n.stage1, n.stage2 = n.stage2,
  #     N.dose = N.dose, N.arm = N.arm, eff.prob = eff.alt,
  #     eff.lower = eff.lower, Cs2 = Cs2, cutoff.step = cutoff.step,
  #     Ce.1 = Ce.1, Ce.2 = Ce.2, type = 2, find.cutoff = F, upfront = upfront
  #   )
  #   power.min <- run.alt.min["Power"]
  #   SR.min <- run.alt.min["SR"]
  # 
  #   tol.val <- if (step == 1) c(0.05, 0.02) else c(0.10, 0.05)
  # 
  #   if (power.min >= target.Power.min - tol.val[1] / target.Power.min &
  #     SR.min >= target.SR.min - tol.val[1] / target.SR.min) {
  #     if (power.min >= target.Power.min - tol.val[2] / target.Power.min &
  #       SR.min >= target.SR.min - tol.val[2] / target.SR.min) {
  #       lower_idx <- rough.idx[val - 1]
  #     } else {
  #       lower_idx <- i
  #     }
  #     break()
  #   }
  # }
  lower_idx <- 1
  if (length(lower_idx) == 0) lower_idx <- 1
  Ce.1 <- Ce.2 <- power.min <- SR.min <- c()
  for (i in lower_idx:dim(sample.sizes)[1]) {
    if (i == dim(sample.sizes)[1]) {
      opt.idx <- FALSE
    }
    n.stage1 <- sample.sizes$n.stage1[i]
    n.stage2 <- sample.sizes$n.stage2[i]

    run.null.max <- one.search(
      n.stage1 = n.stage1, n.stage2 = n.stage2,
      N.dose = N.dose, N.arm = N.arm, eff.prob = eff.null,
      eff.lower = eff.lower, Cs2 = Cs2, cutoff.step = cutoff.step,
      typeI.level = typeI.level, NSR.level = NSR.level, type = 1, upfront = upfront
    )
    Ce.1[i] <- run.null.max["Ce.1"]
    Ce.2[i] <- run.null.max["Ce.2"]

    # run.alt.max <- one.search(
    #   fda.case = fda.case, n.stage1 = n.stage1, n.stage2 = n.stage2,
    #   N.dose = N.dose, N.arm = N.arm, eff.prob = eff.alt,
    #   Ce.1 = Ce.1[i], Ce.2 = Ce.2[i], type = 1, find.cutoff = F, upfront = upfront
    # )
    # power.max[i] <- run.alt.max["Power"]
    # SR.max[i] <- run.alt.max["SR"]

    ### min eff
    run.alt.min <- one.search(
      n.stage1 = n.stage1, n.stage2 = n.stage2,
      N.dose = N.dose, N.arm = N.arm, eff.prob = eff.alt,
      eff.lower = eff.lower, Cs2 = Cs2, cutoff.step = cutoff.step,
      Ce.1 = Ce.1[i], Ce.2 = Ce.2[i], type = 2, find.cutoff = F, upfront = upfront
    )

    power.min[i] <- run.alt.min["Power"]
    SR.min[i] <- run.alt.min["SR"]

    if (power.min[i] >= target.Power.min &
      SR.min[i] >= target.SR.min) {
      break()
    }
  }

  sample.size.table <- data.frame(
    Ce1 = Ce.1, Ce2 = Ce.2,
    Power.min = round(power.min, 4), SR.min = round(SR.min, 4)
  )
  # sample.size.table <- round(sample.size.table, digits = 4, )
  sample.size.table <- cbind(sample.sizes[1:dim(sample.size.table)[1], ], sample.size.table)
  return(list(
    opt.idx = opt.idx,
    sample.size.table = sample.size.table,
    optimal = sample.size.table[nrow(sample.size.table), ]
  ))
}

# # Searching =================================
#
# start <- Sys.time()
# calibration.res <- all.search(
#   fda.case = 1, n.stage1.lower = 10, n.stage1.upper = 40,
#   n.stage2.lower = 2, n.stage2.upper = 40, step = 1,
#   N.dose = 2, N.arm = 4, typeI.level = 0.10, NSR.level = 0.05,
#   target.Power.min = 0.80, target.SR.min = 0.65,
#   eff.lower = 0.25, Cs2 = 0.80, cutoff.step = 1,
#   diff.val = Inf, upfront = TRUE
# )
# end <- Sys.time()
# end - start
#
# calibration.res$optimal

# write to google sheet
# library(googlesheets4)
# google.tab <- gs4_create("sample.size.table", sheets = c("sample.size"))
# sheet_write(data = sample.size.table, ss = google.tab, sheet = "sample.size")

# gs4_find("sample.size.table") %>%
#   googledrive::drive_trash()
