# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Most of this code is from Moka et al. (2024), any modifications   #
# start and end with                                                # 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         START: HUANG ET AL.                       #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                                                                   #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         END: HUANG ET AL.                         #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
library(future)      
library(furrr)       
library(MASS)
library(dplyr)       
library(purrr)    
library(spam)
library(MASS)
library(robustbase)
library(FRBmodelselection)
library(cPCG)


noise.compute.from.SNR = function(SNR,beta,sigma) {
  sigma.2 = t(matrix(beta,ncol=1))%*%sigma%*%matrix(beta,ncol=1)/SNR
  return(sigma.2)
}

cov.X = function(p,c){
  mat.cov = diag(p)
  for (i in 1:p){
    for (j in i:p){
      mat.cov[i,j] = c**(abs(i-j))
      mat.cov[j,i] = mat.cov[i,j]
    }
  }
  return(mat.cov)
}


map.w.to.t = function(w) 1-exp(-w**2)
map.t.to.w = function(t) (-log(1-t))**(0.5)


map_to_combss_object = function(x,tau=0.5,n,X,y){
  t = x[[1]]
  s = as.vector((t>tau))
  beta.hat.s = rep(0,dim(X)[1])
  
  if(sum(s)==0|sum(s)>n){beta.hat.s = 0}else{
    beta.hat.s = solve((t(X[,s])%*%X[,s]))%*%t(X[,s])%*%matrix(y,ncol=1)}
  
  index.cov = which(s==TRUE)
  result = list(index.cov=index.cov,s=s,t=t,Niter=x[[4]],beta.hat.s=beta.hat.s,convergence=x[[3]])
  class(result) = "COMBSS"
  return(result)
}



predict.COMBSS = function(object,newdata){
  if(sum(object$s)==0|sum(object$beta.hat.s)==0){y.pred = rep(0,dim(newdata)[1])}else{
    Xpred = newdata[,object$s] 
    y.pred = Xpred%*%object$beta.hat.s}
  return(y.pred)
}


ADAM.COMBSS = function(X,y,delta=n,lambda,tau=0.5,Niter=200,alpha=0.001,
                       psy=c(0.9,0.999),epoch=10,tol=0.0001,CG=TRUE,trunc=0.001,
                       subsetmap = "V2",
                       criterion = c("MSE", "FAMM"),
                       MMweights = F){
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                         START: HUANG ET AL.                       #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  if (criterion == "FAMM" | criterion == "FRB" ){
    R = rlm(y ~ X, method = "MM", maxit = 500)
    W = diag(R$w)
    X = W^{1/2} %*% X # Rescale the data, based on the weights from the MM-estimator.
    y = W^{1/2} %*% y
  }else{
    R = data.frame(w = NULL, residuals = NULL)
  }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                         END: HUANG ET AL.                         #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  
  
  n = dim(X)[1]
  p = dim(X)[2]
  indice = 1:p
  Xnew = X[,indice]
  XtXnew  = (1/n)*crossprod(Xnew)
  Xtynew = (1/n)*crossprod(Xnew,matrix(y,ncol=1))
  
  t0 = rep(0.5,p)
  w = map.t.to.w(t0)
  c = 10e-8
  #w = w - alpha*u/sqrt(v+c)
  u = rep(0,p)
  v = rep(0,p)
  u.tilde = rep(0,p)
  v.tilde = rep(0,p)
  t.new = t0
  i = 0
  counter = 0
  pnew = p
  wnew = w
  p.new.iter = p
  
  # ADD TPATH
  tpath = vector("list", Niter)
  
  
  while (i<Niter&counter< epoch){
    
    if(pnew<p.new.iter){
      Xnew = X[,indice]
      pnew = length(indice)
      XtXnew  = (1/n)*crossprod(Xnew)
      Xtynew = (1/n)*crossprod(Xnew,matrix(y,ncol=1))
    }
    wnew = w[indice]
    p.new.iter = pnew
    i = i+1
    t.old = t.new 
    gradw = gradient.Combss.Im(w=wnew,p=pnew,X=Xnew,XtXn=XtXnew,Xtyn=Xtynew,delta=n,lambda=lambda,CG=CG,n=n)
    
    u[indice] =psy[1]*u[indice]-(1-psy[1])*gradw
    v[indice] =psy[2]*v[indice]+(1-psy[2])*(gradw*gradw)
    
    u.tilde[indice] = u[indice]/(1-psy[1]**i)
    v.tilde[indice] = v[indice]/(1-psy[2]**i)
    w[indice] = w[indice] + alpha*u.tilde[indice]/sqrt(v.tilde[indice]+c)
    
    t.new = map.w.to.t(w)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #                         START: HUANG ET AL.                       #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    tpath[[i]] = t.new # TPATH
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #                         END: HUANG ET AL.                         #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    
    
    if(!is.null(trunc)){
      if(any(t.new[indice]<trunc)){
        indice = which(t.new>trunc)
        t.new[-c(indice)] = 0
        pnew = length(indice)
      }
      if(is.null(pnew)) break
    }
    
    maxnorm = max(abs(t.old-t.new))
    if(maxnorm < tol){counter = counter + 1} else counter =0
  }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                         START: HUANG ET AL.                       #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  tpath = tpath[1:i]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                         END: HUANG ET AL.                         #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  
  if(i<Niter) convergence = TRUE else convergence = FALSE
  
  if (subsetmap == "V1"){
    t = map.w.to.t(w)
    s = as.vector((t>tau))
    if(sum(s)==0|sum(s)>n) beta.hat.s = 0 else{
      beta.hat.s = solve((t(X[,s])%*%X[,s]))%*%t(X[,s])%*%matrix(y,ncol=1)}  
    
    result = list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,
                  beta.hat.s=beta.hat.s,convergence=convergence,
                  tpath      = tpath[seq_len(i)])
    class(result) = "COMBSS"
    return(result)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    #                         START: HUANG ET AL.                       #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  }else if (subsetmap == "V2"){
    M = SubsetMapV2(tpath = tpath, X = X, y = y, criterion = criterion,
                    weights = R$w, residuals = R$residuals,
                    MMweights = MMweights)  
    return(M)
  }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                         END: HUANG ET AL.                         #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         START: HUANG ET AL.                       #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
SubsetMapV2 = function(X,y,tpath, q = min(nrow(X), ncol(X)) - 1,
                       criterion = c("FAMM", "MSE", "FRB"),
                       weights = R$w, residuals = R$residuals,
                       MMweights){
  # Solve mSE
  critMSE = function(s){
    s = s == 1
    Xs = X[,s]
    beta = solve(t(Xs) %*% Xs) %*% t(Xs) %*% y
    mean((y - Xs %*% beta)^2)
  }
  
  critFAMM = function(s, residuals, weights) {
    critFamm(X      = X, y = y,
             A      = matrix(s, nrow = 1),
             strata = 8,
             B      = 100,
             m      = nrow(X) / 2,
             params = "SM",
             resid = residuals,
             weights = weights,
             MMweights = MMweights)$MSC5
  }
  
  critFRB = function(s){
    msFRB.C(X = X, y = y, A = matrix(s, nrow = 1),
            B = 100, m = round(nrow(X)/2), 
            params = "SM", crit = T,
            MMweights = MMweights)$MSC.5
  }
  
  p  = length(tpath[[1]])          
  Mk = vector("list", q)            
  for (k in 1:q) {                   
    Mk[[k]] = matrix(0, 0, p)
  }
  
  # tpath is a list of vectors
  # t is one vector
  for (t in tpath) {                
    ord = order(t, decreasing = TRUE)   # indices ranked by t-value
    
    for (k in 1:q) {
      s = integer(p)                
      s[ord[1:k]] = 1               # keep the k largest coordinates
      Mk[[k]] = rbind(Mk[[k]], s)   # add to bucket M_k
    }
  }
  sstark = c()
  
  
  for (k in 1:q){
    Mkq = Mk[[k]] # Proposed models of size k ONLY
    if (unique(Mkq) |> nrow() == 1){ # if all the models compared are the same, there is only one choice - no need to run FAMM.
      bestmodel = 1
      sstark = rbind(sstark, Mkq[bestmodel,])
    }else{
      Mkqq = unique(Mkq)
      if (criterion == "MSE"){
        bestmodel = which.min(apply(Mkqq, 1, function(s){
          critMSE(s)  
        }))
      }else if (criterion == "FAMM"){
        bestmodel = which.min(apply(Mkqq, 1, function(s){
          critFAMM(s, residuals = residuals, weights = weights)
        }))
      }else if (criterion == "FRB"){
        bestmodel = which.min(apply(Mkqq, 1, function(s){
          critFRB(s)
        }))
      }
      sstark = rbind(sstark, Mkqq[bestmodel,])
    }
    # print(k/q)
  }
  return(sstark)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         END: HUANG ET AL.                         #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


gradient.Combss.Im = function(w,p,X,XtXn=XtXn,Xtyn=Xtyn,delta,lambda=2,CG,n){
  t = as.vector(map.w.to.t(w))
  B = t((t%d*%XtXn))
  #Lt = t%d*%B
  
  A = t%d*%B
  if(p==1) Lt=A +(1/n)*(delta*(1-t**2)) else Lt=A +(1/n)*(delta*(diag(1-t**2)))
  
  
  if(CG==TRUE){
    u = t%d*%Xtyn # u 
    if(p<n){
      beta.tilde = cgsolve(Lt,u,maxIter = 100) 
    }else{
      rt=(n/delta)*(1/(1-t**2)) 
      temp0 = rt*(t**2)
      temp0 = rt*t
      B = temp0%d*%t(X)
      
      Lttilde = (1/n)*X%*%(t%d*%B) ### Ltilde_t-I
      diag(Lttilde ) =diag(Lttilde )+1  ### B is  Ltilde_t
      
      temp = X%*%matrix(t*rt*u,ncol=1)
      
      beta.tilde = rt*u-(1/n)*B%*%(cgsolve(Lttilde ,temp,maxIter = 100))
      
    }
  }else{
    Lt.inv = solve(Lt)
    beta.tilde = Lt.inv%*%(t%d*%Xtyn)
  }
  t.beta.tilde = matrix(t*beta.tilde,ncol=1)
  at =  XtXn%*%t.beta.tilde-Xtyn
  bt = at -(delta/n)*(t.beta.tilde)
  
  
  if(CG==TRUE){
    if(p<n){
      ct = cgsolve(Lt,matrix(t*at,ncol=1),maxIter = 100)}else{
        #time = proc.time()
        u = t*at
        rt=(n/delta)*(1/(1-t**2)) 
        temp0 = rt*(t**2)
        temp0 = rt*t
        B = temp0%d*%t(X)
        Lttilde = (1/n)*X%*%(t%d*%B) ### Ltilde_t-I
        diag(Lttilde ) =diag(Lttilde )+1  ### B is  Ltilde_t
        
        temp = X%*%matrix(t*rt*u,ncol=1)
        ct = rt*u-(1/n)*B%*%(cgsolve(Lttilde ,temp,maxIter = 100))
        
      }
  }else{
    ct = Lt.inv%*%(matrix(t*at,ncol=1))
  }
  dt = (XtXn-delta/n*diag(p))%*%(matrix(t*ct,ncol=1))
  zetat = 2*(beta.tilde*(at-dt))-2*(bt*ct)
  gradient.g = (zetat+lambda)*2*w*exp(-w**2)
  return(gradient.g)
}


optLambda = function(X, y, test.split = 0.25, criterion = "FAMM"){
  X = as.matrix(X)
  n = nrow(X)
  # test set index
  idx = sample(c(T,F), test.split*n, replace = T)
  Xtest = X[idx,]
  ytest = y[idx]
  Xtrain = X[-idx,]
  ytrain = y[-idx]
  
  CG = TRUE
  alpha = 0.1
  Niter = 1000
  epoch = 10
  tol = 0.001
  tau=0.5
  epoch = 10
  trunc = 0.001
  lambda.max = sum(y*y)/n
  c.grid = 0.8
  nlambda = 50
  grid.comb = rev(lambda.max*(c.grid^(c(0:(nlambda-1)))))
  mse = rep(0,nlambda)
  FAMM = rep(0,nlambda)
  
  for(j in 1:nlambda){
    lam = as.numeric(grid.comb[j])
    c1 = ADAM.COMBSS(Xtrain,ytrain,delta=dim(Xtrain)[1],lambda=lam,tau=tau,Niter=Niter,alpha=alpha,
                     epoch=epoch,tol=tol,CG=CG,trunc=trunc,
                     subsetmap = "V1", criterion = criterion) #criterion is never used here, actually it is being used to rescale the data. 
    y.pred = as.vector(predict.COMBSS(c1,Xtest))
    if(sum(c1$s)>n){mse[j]= 9999}else{
      mse[j] = mean((ytest-y.pred)**2)}
  }
  
  lambda.min = grid.comb[which.min(mse)]
  return(lambda.min)
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         START: HUANG ET AL.                       #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
checkCorrect = function(A, active, p){
  subsets = 1:(active - 1)
  supersets = (active + 1):(p-1)
  true = c(rep(1, active), rep(0, p - active))
  
  correct = c()
  for (model.idx in 1:(p-1)){
    if (model.idx %in% subsets){
      # If number of predictors is less than active
      # Check if it is a subset of the true model
      correct[model.idx] = all(A[model.idx,] <= true)
    }else if (model.idx %in% supersets){
      # If number of predictors is more than active
      # Check if it is a superset of the true model
      correct[model.idx] = all(true <= A[model.idx,])
    }else if (model.idx == active){
      # If number of predictors is the same as active
      # Check that it is exactly the same
      correct[model.idx] = all(true == A[model.idx,])
    }
  }
  return(correct)
}


run_val = function(X,y){
  CG <- TRUE
  alpha <- 0.1
  Niter <- 1000
  epoch <- 10
  tol <- 0.001
  tau=0.5
  epoch <- 10
  trunc <- 0.001
  
  lambda.max <- sum(y*y)/n
  c.grid <- 0.8
  nlambda <- 50
  grid.comb <- rev(lambda.max*(c.grid^(c(0:(nlambda-1)))))
  
  mse <- rep(0,nlambda)
  nselMSE <- rep(0,nlambda)
  famm <- rep(0,nlambda)
  nselFAMM <- rep(0,nlambda)
  
  tMSE = 0
  tFAMM = 0
  for(j in 1:nlambda){
    lam <- as.numeric(grid.comb[j])
    
    start = Sys.time()
    combssMSE <- ADAM.COMBSS(X,y,delta=dim(X)[1],lambda=lam,tau=tau,Niter=Niter,alpha=alpha,epoch=epoch,tol=tol,CG=CG,trunc=trunc, criterion = "MSE", subsetmap = "V1")
    tMSE = tMSE + (Sys.time() - start)
    
    start = Sys.time()
    combssFAMM <- ADAM.COMBSS(X,y,delta=dim(X)[1],lambda=lam,tau=tau,Niter=Niter,alpha=alpha,epoch=epoch,tol=tol,CG=CG,trunc=trunc, criterion = "FAMM", subsetmap = "V1")
    tFAMM = tFAMM + (Sys.time() - start)
    
    nselMSE[j] <- sum(combssMSE$s)
    nselFAMM[j] <- sum(combssFAMM$s)
    
    y.predMSE <- as.vector(predict.COMBSS(combssMSE,Xtest))
    y.predFAMM <- as.vector(predict.COMBSS(combssFAMM,Xtest))
    
    
    if(sum(combssMSE$s)>n){mse[j]<- 9999}else{
      mse[j] <- mean((ytest-y.predMSE)**2)}
    
    if(sum(combssFAMM$s)>n){famm[j]<- 9999}else{
      famm[j] <- mean((ytest-y.predFAMM)**2)}
  }
  
  
  lambda.minMSE <- grid.comb[which.min(mse)]
  lambda.minFAMM <- grid.comb[which.min(famm)]
  
  
  df <- data.frame(
    grid.comb = grid.comb,
    COMBSS  = mse,  nselMSE  = nselMSE,
    rCOMBSS = famm, nselFAMM = nselFAMM
  ) |>
    pivot_longer(
      cols = c(COMBSS, rCOMBSS),
      names_to  = "method",
      values_to = "value"
    ) |>
    mutate(
      nsel = case_when(
        method == "COMBSS"  ~ nselMSE,
        method == "rCOMBSS" ~ nselFAMM,
      ),
      method = factor(method, levels = c("rCOMBSS", "COMBSS"),
                      labels = c("Robust COMBSS", "COMBSS")))
  
  
  list(tMSE, tFAMM, df = df)
}


plot_val = function(df){
  ggplot(df, aes(x = log(grid.comb), y = value, colour = method)) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(rows = vars(method)) +
    scale_x_continuous(
      name = expression(log(lambda)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(y = "MSE (validation set)") +
    theme_few_grid(base_size = 15) +
    scale_colour_manual(values = cols) +
    theme(axis.ticks.x.top = element_blank(), legend.position = "none")
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                         END: HUANG ET AL.                         #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #