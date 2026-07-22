# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Most of this code is from Moka et al. (2024)                      #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("COMBSS/CFuncs.R")
source("Functions.R")
# Low Dimensional
# https://github.com/benoit-liquet/COMBSS-R-VIGNETTE/blob/main/Low_dimensional_example.md
n <- 100
p<- 20
index.true <- seq(1,p,by=2)
K0 <- 10
betald <- rep(0,p)
betald[index.true] <- 1
True.set <- betald !=0
SNR <- 6
rho <- 0.5
matcovX <- cov.X(p,rho)
sigma.2 <- noise.compute.from.SNR(SNR,betald,matcovX)
noise_var <- as.numeric(sigma.2)
set.seed(140)
Xld <- rmvnorm(n,mean=rep(0,p),sigma=matcovX)
yld <- Xld%*%matrix(betald,ncol=1)+rnorm(n,sd=sqrt(sigma.2))

Ntest <- 5000
Xtest <- rmvnorm(Ntest,mean=rep(0,p),sigma=matcovX)
ytest <- Xtest%*%matrix(betald,ncol=1)+rnorm(n,sd=sqrt(sigma.2))


ld0.1 = contaminate(Xld,yld, 0.1, sigma.mat = matcovX, true_beta = betald)
y0.1ld = ld0.1$y_tilde
X0.1ld = ld0.1 |> dplyr::select(-y_tilde) |> as.matrix()

ld0.2 = contaminate(Xld,yld, 0.2, sigma.mat = matcovX, true_beta = betald)  
y0.2ld = ld0.2$y_tilde
X0.2ld = ld0.2 |> dplyr::select(-y_tilde) |> as.matrix()

rld0 = run_val(Xld,yld)
rld0.1 = run_val(X = X0.1ld,y = y0.1ld)
rld0.2 = run_val(X = X0.2ld,y = y0.2ld)


# Higher Dimensional
n <- 1000
p<- 100
index.true <- seq(1,100,by=10)
K0 <- 10
betahd <- rep(0,p)
betahd[index.true] <- 1
True.set <- betahd !=0
SNR <- 6
rho <- 0.8
matcovX <- cov.X(p,rho)
sigma.2 <- noise.compute.from.SNR(SNR,betahd,matcovX)
noise_var <- as.numeric(sigma.2)
set.seed(140)
Xhd <- rmvnorm(n,mean=rep(0,p),sigma=matcovX)
yhd <- Xhd%*%matrix(betahd,ncol=1)+rnorm(n,sd=sqrt(sigma.2))
Ntest <- 5000
Xtest <- rmvnorm(Ntest,mean=rep(0,p),sigma=matcovX)
ytest <- Xtest%*%matrix(betahd,ncol=1)+rnorm(n,sd=sqrt(sigma.2))

hd0.1 = contaminate(Xhd,yhd, 0.1, sigma.mat = matcovX, true_beta = betahd)
y0.1hd = hd0.1$y_tilde
X0.1hd = hd0.1 |> dplyr::select(-y_tilde) |> as.matrix()

hd0.2 = contaminate(Xhd,yhd, 0.2, sigma.mat = matcovX, true_beta = betahd)  
y0.2hd = hd0.2$y_tilde
X0.2hd = hd0.2 |> dplyr::select(-y_tilde) |> as.matrix()

rhd0 = run_val(Xhd,yhd)
rhd0.1 = run_val(X = X0.1hd,y = y0.1hd)
rhd0.2 = run_val(X = X0.2hd,y = y0.2hd)


save(rld0, rld0.1, rld0.2, rhd0, rhd0.1, rhd0.2, file = "Cb.RData")
