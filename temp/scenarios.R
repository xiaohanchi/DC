

n_by_site <- list(
  # continuous y
  # 4 sites
  bind_cols(ctrl = c(50, 40, 30, 100), trt = c(0, 0, 0, 30)),
  bind_cols(ctrl = c(50, 40, 30, 100)*2, trt = c(0, 0, 0, 30)*2), 
  bind_cols(ctrl = c(50, 40, 30, 100)*5, trt = c(0, 0, 0, 30)*5),
  # 6 sites
  bind_cols(ctrl = c(50, 40, 30, 30, 60, 160), trt = c(0, 0, 0, 0, 0, 40)),
  bind_cols(ctrl = c(50, 40, 30, 30, 60, 160)*2, trt = c(0, 0, 0, 0, 0, 40)*2), 
  bind_cols(ctrl = c(50, 40, 30, 30, 60, 160)*5, trt = c(0, 0, 0, 0, 0, 40)*5),
  
  # binary y
  # 4 sites
  bind_cols(ctrl = c(70, 50, 40, 125), trt = c(0, 0, 0, 40)),
  bind_cols(ctrl = c(70, 50, 40, 125)*2, trt = c(0, 0, 0, 40)*2), 
  bind_cols(ctrl = c(70, 50, 40, 125)*5, trt = c(0, 0, 0, 40)*5),
  # 6 sites
  bind_cols(ctrl = c(60, 50, 40, 40, 80, 240), trt = c(0, 0, 0, 0, 0, 60)),
  bind_cols(ctrl = c(60, 50, 40, 40, 80, 240)*2, trt = c(0, 0, 0, 0, 0, 60)*2), 
  bind_cols(ctrl = c(60, 50, 40, 40, 80, 240)*5, trt = c(0, 0, 0, 0, 0, 60)*5)
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
  # 4 sites w/ 5 periods
  rep(0, 5),
  seq(-0.6, 0, length.out = 5),
  c(0, 0.3, 0.6, 0.3, 0),
  # 6 sites w/ 8 periods
  rep(0, 8),
  seq(-0.7, 0, length.out = 8),
  c(0, 0.2, 0.5, 0.8, 0.6, 0.4, 0.2, 0),
  # binary y
  # 4 sites w/ 5 periods
  rep(0, 5),
  seq(-0.8, 0, length.out = 5),
  c(0, 0.5, 1, 0.5, 0),
  # 6 sites w/ 8 periods
  rep(0, 8),
  seq(-1.05, 0, length.out = 8),
  c(0, 0.3, 0.6, 1, 0.75, 0.5, 0.25, 0)
)


active_time <- list(
  # 4 sites w/ 5 periods
  list(
    site1 = 1:3,
    site2 = 2:4,
    site3 = 4:5,
    site4_ctrl = 1:5,
    site4_trt = 5
  ), 
  # 6 sites w/ 8 periods
  list(
    site1 = 1:3,
    site2 = 2:4,
    site3 = 4:5,
    site4 = 5:6,
    site5 = 5:8,
    site6_ctrl = 1:8,
    site6_trt = 8
  )
)


