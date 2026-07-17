# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("COMBSS/CFuncs.R")
source("Functions.R")
combss = function(i, n, p, active,
                  snr = 2, tau,
                  klev, kslo, MMweights = F) {
  Rcpp::sourceCpp("CppFuncs.cpp")
  set.seed(i)
  Niter    = 200        
  alpha    = 0.001      
  epoch    = 10         
  tol      = 1e-4       
  CG       = TRUE       
  trunc    = 0.001   
  data1 = generateData(n = n, p = p, snr = snr, active = active,
                       MSB = FALSE, tau = tau, klev = klev, kslo = kslo)
  
  X = data1$data |> dplyr::select(-y) |> as.matrix()
  y = data1$data$y
  
  lfamm = optLambda(X = X, y = y, criterion = "FAMM")
  
  # optLambda doesn't use criterion specifically, only for FAMM to rescale data
  lmse = lfrb = optLambda(X = X, y = y, criterion = "MSE")
  
  start = Sys.time()
  Afamm1 = ADAM.COMBSS(
    X, y, delta = nrow(X), lambda = lfamm,
    Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
    CG = CG, trunc = trunc, criterion = "FAMM", MMweights = MMweights,
  )
  TimeFamm = as.numeric(Sys.time() - start, units = "secs")
  
  start = Sys.time()
  Amse1 = ADAM.COMBSS(
    X, y, delta = nrow(X), lambda = lmse, 
    Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
    CG = CG, trunc = trunc, criterion = "MSE"
  )
  TimeMSE = as.numeric(Sys.time() - start, units = "secs")
  
  start = Sys.time()
  Afrb1 = ADAM.COMBSS(
    X, y, delta = nrow(X), lambda = lfrb,
    Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
    CG = CG, trunc = trunc, criterion = "FRB", MMweights = MMweights,
  )
  TimeFRB = as.numeric(Sys.time() - start, units = "secs")
  
  
  list(                               # return both logical vectors
    famm = checkCorrect(Afamm1, active, p),
    mse = checkCorrect(Amse1, active, p),
    frb = checkCorrect(Afrb1, active, p),
    TimeFAMM = TimeFamm, TimeMSE  = TimeMSE, TimeFRB = TimeFRB
  )
}

repeats = 1000
plan(multisession, workers = 50)
r1.0 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 100, klev = 2, tau = 0,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm1.0 = r1.0 |> map("famm") |> reduce(`+`)
mse1.0 = r1.0 |> map("mse") |> reduce(`+`)
frb1.0 = r1.0 |> map("frb") |> reduce(`+`)

famm1.0t = r1.0 |> map("TimeFAMM") |> unlist()
mse1.0t = r1.0 |> map("TimeMSE") |> unlist()
frb1.0t = r1.0 |> map("TimeFRB") |> unlist()

r1.0.1 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 100, klev = 2, tau = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm1.0.1 = r1.0.1 |> map("famm") |> reduce(`+`)
mse1.0.1 = r1.0.1 |> map("mse") |> reduce(`+`)
frb1.0.1 = r1.0.1 |> map("frb") |> reduce(`+`)

famm1.0.1t = r1.0.1 |> map("TimeFAMM") |> unlist()
mse1.0.1t = r1.0.1 |> map("TimeMSE") |> unlist()
frb1.0.1t = r1.0.1 |> map("TimeFRB") |> unlist()

save(famm1.0,    mse1.0,    frb1.0,
     famm1.0.1,  mse1.0.1,  frb1.0.1,
     famm1.0t,    mse1.0t,    frb1.0t,
     famm1.0.1t,  mse1.0.1t,  frb1.0.1t,
     file = "COMBSS1.RData")

r2.0 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 100, klev = 2, tau = 0,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.0 = r2.0 |> map("famm") |> reduce(`+`)
mse2.0 = r2.0 |> map("mse") |> reduce(`+`)
frb2.0 = r2.0 |> map("frb") |> reduce(`+`)

famm2.0t = r2.0 |> map("TimeFAMM") |> unlist()
mse2.0t = r2.0 |> map("TimeMSE") |> unlist()
frb2.0t = r2.0 |> map("TimeFRB") |> unlist()

r2.0.1 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 100, klev = 2, tau = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.0.1 = r2.0.1 |> map("famm") |> reduce(`+`)
mse2.0.1 = r2.0.1 |> map("mse") |> reduce(`+`)
frb2.0.1 = r2.0.1 |> map("frb") |> reduce(`+`)

famm2.0.1t = r2.0.1 |> map("TimeFAMM") |> unlist()
mse2.0.1t = r2.0.1 |> map("TimeMSE") |> unlist()
frb2.0.1t = r2.0.1 |> map("TimeFRB") |> unlist()


save(famm2.0,    mse2.0,    frb2.0,
     famm2.0.1,  mse2.0.1,  frb2.0.1,
     famm2.0t,   mse2.0t,   frb2.0t,
     famm2.0.1t, mse2.0.1t, frb2.0.1t,
     file = "COMBSS2.RData")
plan(sequential)
