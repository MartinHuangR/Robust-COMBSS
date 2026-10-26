# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.1 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("Functions.R")
runMS = function(i, n, p, active,
                 snr = 2, tau, m = n/2, B = 100,
                 shift, spread, method, cov = "toeplitz", strata = 8,
                 params = "SM", MMweights = F, base_seed = 1) {
  Rcpp::sourceCpp("CppFuncs.cpp")
  set.seed(i)
  d = generateData(n = n, p = p, active = active, snr = snr,
                   cov = cov, tau = 0, klev = 0, kslo = 0, MSB = F)
  d1 = d$data
  X = d1 |> dplyr::select(-y) |> as.matrix()
  y = d1$y
  beta = d$params$beta
  signal = as.numeric(X %*% beta)
  sigma = sqrt(mean(signal^2)/snr)
  contaminated = sample(1:n, size = floor(tau * n), replace = F)
  if (length(contaminated) > 0){
    y[contaminated] = signal[contaminated] + sigma * (
      shift * sample(c(1), length(contaminated), replace = T) +
        rnorm(length(contaminated), sd = spread))
  }
  A = UOcandModels(X = X, y = y, true_beta = beta, num_candmodels = 5)
  if (!is.numeric(A)){return("Impossible Candidates")}
  switch(method,
         "FAMM" = msFAMMcpp(X = X,
                            y = y,
                            A = A,
                            m = m,
                            strata = strata,
                            B = B,
                            params = params,
                            MMweights = MMweights),
         "FRB" = msFRB.C(X = X,
                         y = y,
                         A = A,
                         m = m,
                         B = B,
                         params = params,
                         MMweights = MMweights),
         "MM" = msMMcpp(X = X,
                        y = y,
                        A = A,
                        m = m,
                        strata = strata,
                        B = B,
                        params = params,
                        MMweights = MMweights))
}

repeats = 1000
# Adjust workers to suit the machine; use repeats = 5 for a pilot run.
plan(multisession, workers = 50)   
configs1 = list(
  list(sim = "Clean", n = 500, p = 50, active = 3, tau = 0, shift = 0, spread = 1, method = "FAMM"),
  list(sim = "Clean", n = 500, p = 50, active = 3, tau = 0, shift = 0, spread = 1, method = "FRB"),
  
  list(sim = "S1", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 1, method = "FAMM"),
  list(sim = "S1", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 1, method = "FRB"),
  list(sim = "S2", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 1, method = "FAMM"),
  list(sim = "S2", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 1, method = "FRB"),
  
  list(sim = "S3", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 0.1, method = "FAMM"),
  list(sim = "S3", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 0.1, method = "FRB"),
  list(sim = "S4", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 0.1, method = "FAMM"),
  list(sim = "S4", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 0.1, method = "FRB"),
  
  
  list(sim = "S5", n = 500, p = 50, active = 3, tau = 0.3, shift = 6, spread = 0.1, method = "FAMM"),
  list(sim = "S5", n = 500, p = 50, active = 3, tau = 0.3, shift = 6, spread = 0.1, method = "FRB")
)

configs2 = list(
  list(sim = "Clean", n = 1000, p = 100, active = 10, tau = 0, shift = 0, spread = 1, method = "FAMM"),
  list(sim = "Clean", n = 1000, p = 100, active = 10, tau = 0, shift = 0, spread = 1, method = "FRB"),
  
  list(sim = "S1", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 1, method = "FAMM"),
  list(sim = "S1", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 1, method = "FRB"),
  list(sim = "S2", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 1, method = "FAMM"),
  list(sim = "S2", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 1, method = "FRB"),
  
  list(sim = "S3", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 0.1, method = "FAMM"),
  list(sim = "S3", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 0.1, method = "FRB"),
  list(sim = "S4", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 0.1, method = "FAMM"),
  list(sim = "S4", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 0.1, method = "FRB"),
  
  
  list(sim = "S5", n = 1000, p = 100, active = 10, tau = 0.3, shift = 6, spread = 0.1, method = "FAMM"),
  list(sim = "S5", n = 1000, p = 100, active = 10, tau = 0.3, shift = 6, spread = 0.1, method = "FRB")
)

configs3 = list(
  list(sim = "Clean", n = 1000, p = 100, active = 10, tau = 0, shift = 0, spread = 1, method = "MM"),
  list(sim = "Clean", n = 500, p = 50, active = 3, tau = 0, shift = 0, spread = 1, method = "MM"),
  
  list(sim = "S1", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 1, method = "MM"),
  list(sim = "S1", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 1, method = "MM"),
  list(sim = "S2", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 1, method = "MM"),
  list(sim = "S2", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 1, method = "MM"),
  
  list(sim = "S3", n = 1000, p = 100, active = 10, tau = 0.1, shift = 4, spread = 0.1, method = "MM"),
  list(sim = "S3", n = 500, p = 50, active = 3, tau = 0.1, shift = 4, spread = 0.1, method = "MM"),
  list(sim = "S4", n = 1000, p = 100, active = 10, tau = 0.2, shift = 4, spread = 0.1, method = "MM"),
  list(sim = "S4", n = 500, p = 50, active = 3, tau = 0.2, shift = 4, spread = 0.1, method = "MM"),
  
  
  list(sim = "S5", n = 1000, p = 100, active = 10, tau = 0.3, shift = 6, spread = 0.1, method = "MM"),
  list(sim = "S5", n = 500, p = 50, active = 3, tau = 0.3, shift = 6, spread = 0.1, method = "MM")
)
configs = c(configs1,configs2,configs3)
results = list()

for (j in seq_along(configs)){
  total = length(configs)
  cfg = configs[[j]]
  message(sprintf("[%d/%d] Running...", j, total))
  r1 = future_map(
    1:repeats,
    runMS,
    n = cfg$n, p = cfg$p, active = cfg$active,
    tau = cfg$tau, shift = cfg$shift, spread = cfg$spread,
    cov = "toeplitz", strata = 8,
    method = cfg$method,
    MMweights = T,
    .options  = furrr_options(seed = TRUE),           
    .progress = TRUE                                  
  ) |> clean()
  r1 = r1 |> mutate(Tau = cfg$tau,
                    n = cfg$n,
                    p = cfg$p,
                    active = cfg$active,
                    shift = cfg$shift,
                    spread = cfg$spread,
                    method = cfg$method,
                    Scenario = cfg$sim)
  new_name = paste0(cfg$method,"t" , cfg$tau, cfg$sim, "n", cfg$n, "UO")
  assign(new_name, r1, envir = .GlobalEnv)
  results[[j]] = r1
  message(sprintf("  -> saved as '%s'", new_name))
}
plan(sequential)
summary_table = bind_rows(results)

save(FAMMt0Cleann500UO, FRBt0Cleann500UO,
     FAMMt0.1S1n500UO, FRBt0.1S1n500UO,
     FAMMt0.2S2n500UO, FRBt0.2S2n500UO,
     FAMMt0.1S3n500UO, FRBt0.1S3n500UO,
     FAMMt0.2S4n500UO, FRBt0.2S4n500UO,
     FAMMt0.3S5n500UO, FRBt0.3S5n500UO,
     
     FAMMt0Cleann1000UO, FRBt0Cleann1000UO,
     FAMMt0.1S1n1000UO, FRBt0.1S1n1000UO,
     FAMMt0.2S2n1000UO, FRBt0.2S2n1000UO,
     FAMMt0.1S3n1000UO, FRBt0.1S3n1000UO,
     FAMMt0.2S4n1000UO, FRBt0.2S4n1000UO,
     FAMMt0.3S5n1000UO, FRBt0.3S5n1000UO,
     
     MMt0Cleann500UO,
     MMt0.1S1n500UO,
     MMt0.2S2n500UO,
     MMt0.1S3n500UO,
     MMt0.2S4n500UO,
     MMt0.3S5n500UO,
     
     MMt0Cleann1000UO,
     MMt0.1S1n1000UO,
     MMt0.2S2n1000UO,
     MMt0.1S3n1000UO,
     MMt0.2S4n1000UO,
     MMt0.3S5n1000UO,
     summary_table, file = "RMS/MSC-NS.RData")

