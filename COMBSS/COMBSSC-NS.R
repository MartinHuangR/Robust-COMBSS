# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("CFuncs.R")
source("Functions.R")
combss = function(i, n, p, active,
                  snr = 2, tau, shift, spread,
                  klev, kslo, MMweights = F) {
  Rcpp::sourceCpp("CppFuncs.cpp")
  set.seed(i)
  Niter    = 200        
  alpha    = 0.001      
  epoch    = 10         
  tol      = 1e-4       
  CG       = TRUE       
  trunc    = 0.001   
  # X is unchanged; only y is contaminated (as in the MS simulation)
  data1 = generateData(n = n, p = p, snr = snr, active = active, cov = "toeplitz",
                       MSB = FALSE, tau = 0, klev = klev, kslo = kslo)
  
  X = data1$data |> dplyr::select(-y) |> as.matrix()
  y = data1$data$y
  beta = data1$params$beta
  signal = as.numeric(X %*% beta)
  sigma = sqrt(mean(signal^2)/snr)
  contaminated = sample(1:n, size = floor(tau * n), replace = F)
  if (length(contaminated) > 0){
    y[contaminated] = signal[contaminated] + sigma * (
      shift * sample(c(1), length(contaminated), replace = T) + # non symmetric error
        rnorm(length(contaminated), sd = spread))
  }
  
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

#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#
# Configuration 1: n = 200, p = 20, active = 3
#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#

# Clean
r1.Clean = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0, shift = 0, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.Clean = r1.Clean |> map("famm") |> reduce(`+`)
mse1.Clean = r1.Clean |> map("mse") |> reduce(`+`)
frb1.Clean = r1.Clean |> map("frb") |> reduce(`+`)

famm1.Cleant = r1.Clean |> map("TimeFAMM") |> unlist()
mse1.Cleant = r1.Clean |> map("TimeMSE") |> unlist()
frb1.Cleant = r1.Clean |> map("TimeFRB") |> unlist()

r1.S1 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0.1, shift = 4, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.S1 = r1.S1 |> map("famm") |> reduce(`+`)
mse1.S1 = r1.S1 |> map("mse") |> reduce(`+`)
frb1.S1 = r1.S1 |> map("frb") |> reduce(`+`)

famm1.S1t = r1.S1 |> map("TimeFAMM") |> unlist()
mse1.S1t = r1.S1 |> map("TimeMSE") |> unlist()
frb1.S1t = r1.S1 |> map("TimeFRB") |> unlist()

r1.S2 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0.2, shift = 4, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.S2 = r1.S2 |> map("famm") |> reduce(`+`)
mse1.S2 = r1.S2 |> map("mse") |> reduce(`+`)
frb1.S2 = r1.S2 |> map("frb") |> reduce(`+`)

famm1.S2t = r1.S2 |> map("TimeFAMM") |> unlist()
mse1.S2t = r1.S2 |> map("TimeMSE") |> unlist()
frb1.S2t = r1.S2 |> map("TimeFRB") |> unlist()

r1.S3 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0.1, shift = 4, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.S3 = r1.S3 |> map("famm") |> reduce(`+`)
mse1.S3 = r1.S3 |> map("mse") |> reduce(`+`)
frb1.S3 = r1.S3 |> map("frb") |> reduce(`+`)

famm1.S3t = r1.S3 |> map("TimeFAMM") |> unlist()
mse1.S3t = r1.S3 |> map("TimeMSE") |> unlist()
frb1.S3t = r1.S3 |> map("TimeFRB") |> unlist()

r1.S4 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0.2, shift = 4, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.S4 = r1.S4 |> map("famm") |> reduce(`+`)
mse1.S4 = r1.S4 |> map("mse") |> reduce(`+`)
frb1.S4 = r1.S4 |> map("frb") |> reduce(`+`)

famm1.S4t = r1.S4 |> map("TimeFAMM") |> unlist()
mse1.S4t = r1.S4 |> map("TimeMSE") |> unlist()
frb1.S4t = r1.S4 |> map("TimeFRB") |> unlist()

r1.S5 = future_map(
  1:repeats,
  combss,
  n = 200, p = 20, active = 3,
  kslo = 0, klev = 0, tau = 0.3, shift = 6, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
famm1.S5 = r1.S5 |> map("famm") |> reduce(`+`)
mse1.S5 = r1.S5 |> map("mse") |> reduce(`+`)
frb1.S5 = r1.S5 |> map("frb") |> reduce(`+`)

famm1.S5t = r1.S5 |> map("TimeFAMM") |> unlist()
mse1.S5t = r1.S5 |> map("TimeMSE") |> unlist()
frb1.S5t = r1.S5 |> map("TimeFRB") |> unlist()

save(famm1.Clean,           mse1.Clean,           frb1.Clean,
     famm1.S1,       mse1.S1,       frb1.S1,
     famm1.S2,       mse1.S2,       frb1.S2,
     famm1.S3,      mse1.S3,      frb1.S3,
     famm1.S4,      mse1.S4,      frb1.S4,
     famm1.S5,  mse1.S5,  frb1.S5,
     famm1.Cleant,          mse1.Cleant,          frb1.Cleant,
     famm1.S1t,      mse1.S1t,      frb1.S1t,
     famm1.S2t,      mse1.S2t,      frb1.S2t,
     famm1.S3t,     mse1.S3t,     frb1.S3t,
     famm1.S4t,     mse1.S4t,     frb1.S4t,
     famm1.S5t, mse1.S5t, frb1.S5t,
     file = "COMBSS1C-NS.RData")

#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#
# Configuration 2: n = 500, p = 50, active = 5
#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#

# Clean
r2.Clean = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0, shift = 0, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.Clean = r2.Clean |> map("famm") |> reduce(`+`)
mse2.Clean = r2.Clean |> map("mse") |> reduce(`+`)
frb2.Clean = r2.Clean |> map("frb") |> reduce(`+`)

famm2.Cleant = r2.Clean |> map("TimeFAMM") |> unlist()
mse2.Cleant = r2.Clean |> map("TimeMSE") |> unlist()
frb2.Cleant = r2.Clean |> map("TimeFRB") |> unlist()

r2.S1 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0.1, shift = 4, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.S1 = r2.S1 |> map("famm") |> reduce(`+`)
mse2.S1 = r2.S1 |> map("mse") |> reduce(`+`)
frb2.S1 = r2.S1 |> map("frb") |> reduce(`+`)

famm2.S1t = r2.S1 |> map("TimeFAMM") |> unlist()
mse2.S1t = r2.S1 |> map("TimeMSE") |> unlist()
frb2.S1t = r2.S1 |> map("TimeFRB") |> unlist()


r2.S2 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0.2, shift = 4, spread = 1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.S2 = r2.S2 |> map("famm") |> reduce(`+`)
mse2.S2 = r2.S2 |> map("mse") |> reduce(`+`)
frb2.S2 = r2.S2 |> map("frb") |> reduce(`+`)

famm2.S2t = r2.S2 |> map("TimeFAMM") |> unlist()
mse2.S2t = r2.S2 |> map("TimeMSE") |> unlist()
frb2.S2t = r2.S2 |> map("TimeFRB") |> unlist()

r2.S3 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0.1, shift = 4, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.S3 = r2.S3 |> map("famm") |> reduce(`+`)
mse2.S3 = r2.S3 |> map("mse") |> reduce(`+`)
frb2.S3 = r2.S3 |> map("frb") |> reduce(`+`)

famm2.S3t = r2.S3 |> map("TimeFAMM") |> unlist()
mse2.S3t = r2.S3 |> map("TimeMSE") |> unlist()
frb2.S3t = r2.S3 |> map("TimeFRB") |> unlist()

r2.S4 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0.2, shift = 4, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.S4 = r2.S4 |> map("famm") |> reduce(`+`)
mse2.S4 = r2.S4 |> map("mse") |> reduce(`+`)
frb2.S4 = r2.S4 |> map("frb") |> reduce(`+`)

famm2.S4t = r2.S4 |> map("TimeFAMM") |> unlist()
mse2.S4t = r2.S4 |> map("TimeMSE") |> unlist()
frb2.S4t = r2.S4 |> map("TimeFRB") |> unlist()

r2.S5 = future_map(
  1:repeats,
  combss,
  n = 500, p = 50, active = 5,
  kslo = 0, klev = 0, tau = 0.3, shift = 6, spread = 0.1,
  MMweights = T,
  .options  = furrr_options(seed = TRUE),           
  .progress = TRUE                                  
)
famm2.S5 = r2.S5 |> map("famm") |> reduce(`+`)
mse2.S5 = r2.S5 |> map("mse") |> reduce(`+`)
frb2.S5 = r2.S5 |> map("frb") |> reduce(`+`)

famm2.S5t = r2.S5 |> map("TimeFAMM") |> unlist()
mse2.S5t = r2.S5 |> map("TimeMSE") |> unlist()
frb2.S5t = r2.S5 |> map("TimeFRB") |> unlist()

save(famm2.Clean,           mse2.Clean,           frb2.Clean,
     famm2.S1,       mse2.S1,       frb2.S1,
     famm2.S2,       mse2.S2,       frb2.S2,
     famm2.S3,      mse2.S3,      frb2.S3,
     famm2.S4,      mse2.S4,      frb2.S4,
     famm2.S5,  mse2.S5,  frb2.S5,
     famm2.Cleant,          mse2.Cleant,          frb2.Cleant,
     famm2.S1t,      mse2.S1t,      frb2.S1t,
     famm2.S2t,      mse2.S2t,      frb2.S2t,
     famm2.S3t,     mse2.S3t,     frb2.S3t,
     famm2.S4t,     mse2.S4t,     frb2.S4t,
     famm2.S5t, mse2.S5t, frb2.S5t,
     file = "COMBSS2C-NS.RData")
plan(sequential)