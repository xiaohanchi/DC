

n_by_site <- list(
  # continuous y
  # 4 sites
  bind_cols(ctrl = c(30, 40, 50, 60), trt = c(0, 0, 0, 30)),
  bind_cols(ctrl = c(30, 40, 50, 60)*2, trt = c(0, 0, 0, 30)*2), 
  bind_cols(ctrl = c(30, 40, 50, 60)*5, trt = c(0, 0, 0, 30)*5),
  # 6 sites
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100), trt = c(0, 0, 0, 0, 0, 50)),
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100)*2, trt = c(0, 0, 0, 0, 0, 50)*2), 
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100)*5, trt = c(0, 0, 0, 0, 0, 50)*5),
  
  # binary y
  # 4 sites
  bind_cols(ctrl = c(30, 40, 50, 60), trt = c(0, 0, 0, 30)) * 2,
  bind_cols(ctrl = c(30, 40, 50, 60)*2, trt = c(0, 0, 0, 30)*2) * 2, 
  bind_cols(ctrl = c(30, 40, 50, 60)*5, trt = c(0, 0, 0, 30)*5) * 2,
  # 6 sites
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100), trt = c(0, 0, 0, 0, 0, 50)) * 2,
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100)*2, trt = c(0, 0, 0, 0, 0, 50)*2) * 2, 
  bind_cols(ctrl = c(30, 40, 50, 30, 40, 100)*5, trt = c(0, 0, 0, 0, 0, 50)*5) * 2
  
)


site_delta <- list(
  # continuous y
  # 4 sites
  c(0, 0, 0, 0),
  c(-1, 0, 0, 0),
  c(-0.5, 0.5, 0.5, 0),
  # 6 sites
  c(0, 0, 0, 0, 0, 0),
  c(-1, 1, 0, 0, 0, 0),
  c(-0.5, -0.5, 0.5, 0.5, 0, 0),
  
  # binary y
  # 4 sites
  c(0, 0, 0, 0),
  c(-1.2, 0, 0, 0),
  c(-0.8, 0.8, 0.8, 0),
  # 6 sites
  c(0, 0, 0, 0, 0, 0),
  c(-1.2, 1.2, 0, 0, 0, 0),
  c(-0.8, -0.8, 0.8, 0.8, 0, 0)
)

drift_mat <- list(
  # continuous y
  # 4 sites w/ 6 periods
  rep(0, 6),
  seq(-0.6, 0, length.out = 6),
  c(0, 0.3, 0.6, 0.45, 0.15, 0),
  # 6 sites w/ 10 periods
  rep(0, 10),
  seq(-0.81, 0, length.out = 10),
  c(0, 0.2, 0.4, 0.6, 0.8, 0.6, 0.45, 0.3, 0.2, 0),
  
  # binary y
  # 4 sites w/ 6 periods
  rep(0, 6),
  seq(-1, 0, length.out = 6),
  c(0, 0.5, 1, 0.75, 0.25, 0),
  # 6 sites w/ 10 periods
  rep(0, 10),
  seq(-1.35, 0, length.out = 10),
  c(0, 0.25, 0.5, 0.75, 1, 0.8, 0.6, 0.4, 0.2, 0)
)


active_time <- list(
  # 4 sites w/ 6 periods
  list(
    site1 = 1:4,
    site2 = 2:5,
    site3 = 4:6,
    site4_ctrl = 5:6,
    site4_trt = 6
  ), 
  # 6 sites w/ 10 periods
  list(
    site1 = 1:4,
    site2 = 2:5,
    site3 = 4:6,
    site4 = 6:8,
    site5 = 7:10,
    site6_ctrl = 9:10,
    site6_trt = 10
  )
)


