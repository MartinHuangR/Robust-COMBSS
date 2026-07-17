# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# General R-Code functions                                          #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
library(patchwork)
library(ggh4x)
library(furrr)
library(purrr)
library(Rcpp)
library(FRB)
library(FRBmodelselection)
library(pense)
library(MASS)
library(robustbase)
library(dplyr)
library(tidyverse)
contaminate = function(X,y, tau, klev = 2, kslo = 100, sigma.mat, true_beta){
  p = ncol(X)
  n = nrow(X)
  contaminated = sample(1:n, size = floor(tau * n), replace = F)
  
  X_tilde = X
  y_tilde = y
  beta_tilde = true_beta
  beta_tilde[true_beta != 0] = beta_tilde[true_beta != 0] * (1 + kslo)
  beta_tilde[true_beta == 0] = kslo * max(abs(true_beta))
  
  
  for (row in contaminated){
    a_tilde = runif(p,-1,1)
    a = a_tilde - as.numeric(((1/p)*t(a_tilde) %*% rep(1,p)))
    theta = mvnfast::rmvn(n = 1, mu = rep(0,p), sigma = 0.01 * diag(p))
    x_tilde = theta + (klev*a)/as.numeric((sqrt(t(a) %*% solve(sigma.mat) %*% a)))
    X_tilde[row,] = x_tilde
    
    y_tilde[row] = as.numeric(x_tilde %*% beta_tilde)  # Ensures a scalar output
  }
  data.frame(X_tilde, y_tilde)
}


generateData = function(n = NULL, p = NULL, snr = NULL, active = NULL,
                        cov = "toeplitz",
                        tau = NULL, klev = NULL, kslo = NULL,
                        MSB = NULL, error = NULL){
  if (MSB == F){
    error = NULL
    if (cov == "toeplitz"){
      sigma.mat = 0.5^abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))  
      X = mvtnorm::rmvnorm(n  = n, sigma = sigma.mat)
    }else if (cov == "ind"){
      sigma.mat = diag(p)
      X = matrix(rnorm(n * p), nrow = n, ncol = p)
    }
    beta = c(rep(1,active) ,rep(0,p-active))
    signal = sqrt(mean((as.matrix(X) %*% as.matrix(beta))^2))
    sigma = as.numeric(signal/sqrt(snr))
    y = as.matrix(X)%*%as.matrix(beta) + rnorm(nrow(X), 0, sd = sigma)
    if (tau > 0){
      Xy = contaminate(X,y, tau = tau, klev = klev, kslo = kslo, sigma.mat = sigma.mat, true_beta = beta)
      
      X.contam = Xy |> dplyr::select(-y_tilde) |> as.matrix()    
      y.contam = Xy$y_tilde  
    }else{
      X.contam = X
      y.contam = y
    }
  }else if (MSB == T){
    n = 64; p = 2; active = 1; cov = NULL; snr = NULL; error = error
    X1 = runif(64, 0, 1)
    X2 = runif(64, 0 ,1)
    if (error == "3/8"){e = (5/8)*rnorm(64, 0, 1) + (3/8)*rnorm(64, (30 - 2 - 2*X1), 1)}
    if (error == "1/4"){e = (3/4)*rnorm(64, 0, 1) + (1/4)*rnorm(64, (30 - 2 - 2*X1), 1)}
    if (error == "1/8"){e = (7/8)*rnorm(64, 0, 1) + (1/8)*rnorm(64, (30 - 2 - 2*X1), 1)}
    if (error == "N"){e = rnorm(64)}
    if (error == "S"){e = rnorm(64)/runif(64)}
    if (error == "C"){e = rcauchy(64)}
    y.contam = 2 + 2*X1 + 0*X2 + e
    X.contam = cbind(X1, X2)
  }
  out = cbind(as.data.frame(X.contam), y =  y.contam)
  params = list(klev = klev,
                kslo = kslo,
                tau = tau,
                cov = cov,
                active = active,
                snr = snr,
                error = error,
                beta = beta)
  list(data = out,
       params = params)
}

UOcandModels = function(X, y, true_beta, num_candmodels = 5){
  
  # ensure there are enough models.
  underfitting = overfitting = floor((num_candmodels - 1)/2)
  if ((underfitting > sum(true_beta != 0)) ||
      (overfitting + sum(true_beta != 0) >= ncol(X))) {
    return("Too many candidate models")
  }
  tot_cand = underfitting + overfitting + 1
  
  A = matrix(data = 0, nrow = tot_cand, ncol = ncol(X))
  A[1,] = as.numeric(true_beta != 0)
  
  # underfitting
  for (i in 1:underfitting){
    if (i == 1){beta = true_beta}
    
    # find the abs(min) beta and remove it
    non_zero = which(beta != 0)
    max_beta = non_zero[which.min(abs(beta[non_zero]))]# if tie, just take the first located index
    beta[max_beta] = 0 
    A[1+i,] = as.numeric(beta != 0)
  }
  
  # overfitting
  zero_var_idx = which(true_beta == 0)
  corr = abs(cor(X[, zero_var_idx], y))
  ordered = order(corr, decreasing = T) 
  for (i in 1:overfitting){
    if (i == 1){beta = as.numeric(true_beta != 0)}
    max_corr = zero_var_idx[ordered[i]]
    beta[max_corr] = 1
    A[1 + underfitting + i, ] = beta
  }
  A
}

clean = function(r){
  do.call(rbind, r) |> group_by(Criterion) |> summarise(across(c(True, M2, M3, M4, M5), sum),
                                                        TotalTime = sum(Time),
                                                        Time = mean(Time),
                                                        .groups = 'drop')
}


# https://github.com/msalibian/FRBmodelselection/blob/master/R/FRBmodelselection.R
Rho.MM = function(x, cc=4.685061) {
  # Tukey's bisquare rho
  tmp = x/cc
  tmp2 = cc^2*tmp^2*(1-tmp^2+tmp^4/3)/2
  tmp2[ abs(tmp) > 1 ] = cc^2/6
  return(tmp2)
}

msFRB.C = function(X, y, A, B, m, params, crit = F, MMweights = F){
  # This code is from https://github.com/msalibian/FRBmodelselection/tree/master
  start = Sys.time()
  d1 = data.frame(y = y , X)
  xx = model.matrix(y ~., data = d1)
  y = as.vector(y)
  n = nrow(X)
  p = ncol(X)
  
  # Pre-compute bootstrap samples
  boot.samp = matrix(sample(n, B*m, replace = T), B, m)
  if (params == "MSB"){
    # Sigma is the S-scale estimate of the largest model.
    largest_model = apply(A, 1, sum) |> which.max()
    X.largest = A[largest_model,] == 1
    if (!(sum(X.largest) == ncol(X))){
      RLM = MASS::rlm(y ~ as.matrix(X[,X.largest]), method='MM', maxit = 500)  
    }else{
      RLM = MASS::rlm(y ~ as.matrix(X), method='MM', maxit = 500)  
    }
    s.scale = summary(RLM)$sigma # instead of mad.full.ro for sigma. 
  }else{
    RLM = MASS::rlm(y ~ as.matrix(X), method='MM', maxit = 500)  
    # rho function 
    mad.full.ro = mad(RLM$residuals);
    # parameters for robust model selection criterion
    b.huber  = 2; 
  }
  
  if (MMweights){
    if (params == "MSB" && !(sum(X.largest) == ncol(X))){
      RLM.full = MASS::rlm(y ~ as.matrix(X), method='MM', maxit = 500)
      w.crit = as.vector(RLM.full$w)
    }else{
      w.crit = as.vector(RLM$w)
    }
  }else{
    w.crit = rep(1, n)
  }
  
  # Full parameters can be found in FRBmodelselection::rlm.control
  my.ctrl = rlm.control(M = B)
  X.full = cbind(1, X)
  num.models = nrow(A)
  MSC.5 = MSC.2 = c()
  M05.ind = M02.ind = rep(0,num.models)
  A = cbind(1, A) 
  
  
  for (pp in 1:num.models){
    X.red = A[pp,] == 1
    x0 = as.matrix(xx[, X.red], ncol = ncol(xx[, X.red])) # always assume intercept (1st column) 
    tmp = roboot(x=x0, y=y, nboot=m, boot.samp=boot.samp,
                 control=my.ctrl)
    res.boot2 = apply(tmp$ours, 1, function (coefs){
      y - x0 %*% matrix(coefs, ncol = 1)
    })
    
    if (params == "MSB"){
      MSC.2[pp] = apply(res.boot2, 2, function(res){
        sum(w.crit * Rho.MM(res/s.scale))               # (weighted) reduces to sum(Rho.MM(.)) when w_i = 1
      }) |> sum()
    }else{
      MSC.2[pp] = apply(res.boot2, 2, function(res){
        small = as.vector(res)^2 < (b.huber*mad.full.ro)^2  # (added) logical index, same set as original subsetting
        (sum(w.crit[small] * (res[small])^2) + 
            ((b.huber*mad.full.ro)^2)*sum(w.crit[!small]))/(n);  # (weighted) reduces to original when w_i = 1
      }) |> sum()
    }
    pDim     = sum(X.red)
    MSC.5[pp] = MSC.2[pp]/B;
    
    
    res = y - as.matrix(x0) %*% matrix(tmp$coefficients, ncol = 1)
    if (params == "MSB"){
      Ap = sum(X.red) + 1
      penal = (sum(w.crit * Rho.MM(res/s.scale)) +       # (weighted) in-sample loss; penalty term unweighted
                 log(n)*(Ap-1)*(s.scale^2) )/(n) # First part of equation (10)
      MSC.5[pp] = penal + MSC.5[pp];
    }else{
      small = as.vector(res)^2 < (b.huber*mad.full.ro)^2   # (added) logical index, same set as original subsetting
      penal = ((sum(w.crit[small] * (res[small])^2) + 
                  ((b.huber*mad.full.ro)^2)*sum(w.crit[!small]) +
                  2*log(n)*pDim*(mad.full.ro^2))/(n))    # (weighted) loss; penalty unweighted
      MSC.5[pp] = penal + MSC.5[pp];
    }
  }
  end = Sys.time() 
  if (crit == F){
    M05.ind[which.min(MSC.5)] = M05.ind[which.min(MSC.5)] + 1
    M02.ind[which.min(MSC.2)] = M02.ind[which.min(MSC.2)] + 1
    
    r = rbind(M05.ind |> as.data.frame() |> t() |> data.frame(),
              M02.ind |> as.data.frame() |> t() |> data.frame())
    r = tibble::rownames_to_column(r)
    r$Time = end - start
    colnames(r) = c("Criterion", "True", paste0("M",2:num.models), "Time")
    r$Method = "FRB"
    return(r)  
  }else{
    return(data.frame(MSC.2, MSC.5))
  }
}

# Theme function.
theme_few_grid = function (base_size = 12, base_family = "") 
{
  gray <- "#4D4D4D"
  black <- "#000000"
  theme_bw(base_size = base_size, base_family = base_family) + 
    theme(line = element_line(colour = gray), rect = element_rect(fill = "white", 
                                                                  colour = NA), text = element_text(colour = black), 
          axis.ticks = element_line(colour = gray), legend.key = element_rect(colour = NA), 
          panel.border = element_rect(colour = gray), panel.grid = element_line(color = alpha("black", 0.05)), 
          strip.background = element_rect(fill = "#f7f9fc", colour = gray))
}

to_cm <- function(x) {
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]
  ft <- as.numeric(parts[1])
  inch <- as.numeric(parts[2])
  (ft * 12 + inch) * 2.54
}