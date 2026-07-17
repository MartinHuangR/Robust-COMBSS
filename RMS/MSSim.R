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
                 klev, kslo, method, cov = "toeplitz", strata = 8,
                 params = "SM", MMweights = F) {
  Rcpp::sourceCpp("CppFuncs.cpp")
  set.seed(i)
  start = Sys.time()
  d = generateData(n = n, p = p, snr = snr, active = active, cov = cov, tau = tau, klev = klev, kslo = kslo, MSB = F)
  d1 = d$data
  A = UOcandModels(X = d1 |> dplyr::select(-y) |> as.matrix(),
                   true_beta = d$params$beta,
                   y = d1$y)
  if (!is.numeric(A)){return("Impossible Candidates")}
  switch(method,
         "FAMM" = msFAMMcpp(X = d1 |> dplyr::select(-y)|> as.matrix(ncol = p),
                            y = d1$y,
                            A = A,
                            m = m,
                            strata = strata,
                            B = B,
                            params = params,
                            MMweights = MMweights),
         "MM" = msMMcpp(X = d1 |> dplyr::select(-y) |> as.matrix(ncol  = p),
                        y = d1$y,
                        A = A,
                        m = m,
                        strata = strata,
                        B = B,
                        params = params,
                        MMweights = MMweights),
         "FRB" = msFRB.C(X = d1 |> dplyr::select(-y),
                         y = d1$y,
                         A = A,
                         m = m,
                         B = B,
                         params = params,
                         MMweights = MMweights))
}

repeats = 1000
plan(multisession, workers = 60)   
configs1 = list(
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0, kslo = 0, klev = 0, method = "FAMM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0, kslo = 0, klev = 0, method = "FRB"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0, kslo = 0, klev = 0, method = "MM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.05, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.05, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.05, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.1, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.1, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.1, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.15, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.15, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.15, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.2, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.2, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s1", n = 500, p = 50, active = 3, tau = 0.2, kslo = 100, klev = 2, method = "MM")
  
)

for (j in seq_along(configs1)){
  total = length(configs1)
  cfg = configs1[[j]]
  message(sprintf("[%d/%d] Running...", j, total))
  r1 = future_map(
    1:repeats,
    runMS,
    n = cfg$n, p = cfg$p, active = cfg$active,
    kslo = cfg$kslo, klev = cfg$klev, tau = cfg$tau,
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
                    kslo = cfg$kslo,
                    klev = cfg$klev,
                    method = cfg$method)
  new_name = paste0(cfg$method,"t" , cfg$tau, cfg$sim, "UO")
  assign(new_name, r1, envir = .GlobalEnv)
  message(sprintf("  -> saved as '%s'", new_name))
}
plan(sequential)
save(FAMMt0s1UO, FRBt0s1UO, MMt0s1UO,
     FAMMt0.05s1UO, FRBt0.05s1UO, MMt0.05s1UO,
     FAMMt0.1s1UO, FRBt0.1s1UO, MMt0.1s1UO,
     FAMMt0.15s1UO, FRBt0.15s1UO, MMt0.15s1UO,
     FAMMt0.2s1UO, FRBt0.2s1UO, MMt0.2s1UO, file = "MS1.RData")
#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#
repeats = 1000
plan(multisession, workers = 60) 
configs2 = list(
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0, kslo = 0, klev = 0, method = "FAMM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0, kslo = 0, klev = 0, method = "FRB"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0, kslo = 0, klev = 0, method = "MM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.05, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.05, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.05, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.1, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.1, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.1, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.15, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.15, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.15, kslo = 100, klev = 2, method = "MM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.2, kslo = 100, klev = 2, method = "FAMM"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.2, kslo = 100, klev = 2, method = "FRB"),
  list(sim = "s2", n = 1000, p = 100, active = 10, tau = 0.2, kslo = 100, klev = 2, method = "MM")
  
)

for (j in seq_along(configs2)){
  total = length(configs2)
  cfg = configs2[[j]]
  message(sprintf("[%d/%d] Running...", j, total))
  r1 = future_map(
    1:repeats,
    runMS,
    n = cfg$n, p = cfg$p, active = cfg$active,
    kslo = cfg$kslo, klev = cfg$klev, tau = cfg$tau,
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
                    kslo = cfg$kslo,
                    klev = cfg$klev,
                    method = cfg$method)
  new_name = paste0(cfg$method,"t" , cfg$tau, cfg$sim, "UO")
  assign(new_name, r1, envir = .GlobalEnv)
  message(sprintf("  -> saved as '%s'", new_name))
}
plan(sequential)
save(FAMMt0s2UO, FRBt0s2UO, MMt0s2UO,
     FAMMt0.05s2UO, FRBt0.05s2UO, MMt0.05s2UO,
     FAMMt0.1s2UO, FRBt0.1s2UO, MMt0.1s2UO,
     FAMMt0.15s2UO, FRBt0.15s2UO, MMt0.15s2UO,
     FAMMt0.2s2UO, FRBt0.2s2UO, MMt0.2s2UO, file = "MS2.RData")
#--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--#
runMSB = function(errors, repeats = 1000, method, A,
                  strata = NULL, B, m, params = NULL, crit = F,
                  MMweights = F){
  start = Sys.time()
  out = 1:repeats |> map_dfr(function(i){
    set.seed(i + 1)
    d1 = generateData(MSB = TRUE, error = errors)$data
    switch(method,
           "FAMM" = msFAMMcpp(X = d1 |> dplyr::select(-y)|> as.matrix(ncol = p),
                              y = d1$y,
                              A = A,
                              m = m,
                              strata = 8,
                              B = 100,
                              params = params,
                              crit = crit,
                              MMweights = MMweights),
           "MM" = msMMcpp(X = d1 |> dplyr::select(-y) |> as.matrix(ncol = p),
                          y = d1$y,
                          A = A,
                          m = m,
                          strata = 8,
                          params = params,
                          B = 100,
                          MMweights = MMweights),
           "FRB" = msFRB.C(X = d1 |> dplyr::select(-y),
                           y = d1$y,
                           A = A,
                           m = m,
                           B = 100,
                           params = params,
                           crit = crit,
                           MMweights = MMweights))
  },
  .progress = TRUE
  )
  if (crit == F){
    out = out |>
      group_by(Criterion) |>
      summarise(
        across(c(True, M2, M3, M4), sum),
        across(Time, mean),
        .groups = "drop")
    out$TotalTime = Sys.time() - start
    out$Error = errors
    return(out)
  }else{
    out$model = rep(1:4,repeats)
    out$Error = errors
    return(out)
  }
}

A = matrix(rbind(c(1,0), # 1,2 (TRUE)
                 c(0,0), # 1
                 c(0,1), # 1,3
                 c(1,1)),# 1,2,3
           ncol = 2, nrow = 4)
repeats = 1000
errors = c("N", "1/8", "1/4", "3/8", "S", "C")
library(furrr)
library(purrr)
library(future)

FAMM = map_dfr(
  errors,
  runMSB,
  repeats = 1000,
  method = "FAMM",
  A      = A,
  m      = 24,
  params = "SM",
  MMweights = T,
  .progress = TRUE                         
)

MM = map_dfr(
  errors,
  runMSB,
  repeats = 1000,
  method = "MM",
  A      = A,
  m      = 24,
  params = "SM",
  MMweights = T,
  .progress = TRUE                         
)


FRB = map_dfr(
  errors,
  runMSB,
  repeats = 1000,
  method = "FRB",
  A      = A,
  m      = 24,
  params = "SM",
  MMweights = T,
  .progress = TRUE                         
)

# load("MSsm.RData")
save(FAMM, FRB, MM, file = "MSsm.RData")
