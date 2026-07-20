# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.3 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("COMBSS/CFuncs.R")
source("Functions.R")
Rcpp::sourceCpp("CppFuncs.cpp")
dq = read.csv("NBA/NBAOutlier.csv")

dq = dq[,2:ncol(dq)]
y = dq$ast
X = dq  |> dplyr::select(-c(ast)) |> scale() |> as.matrix() 
set.seed(1) 

Niter    = 200        
alpha    = 0.001      
epoch    = 10         
tol      = 1e-4       
CG       = TRUE       
trunc    = 0.001   

lfamm = optLambda(X = X, y = y, criterion = "FAMM")
lmse = lfrb = optLambda(X = X, y = y, criterion = "MSE")
start = Sys.time()
Afamm = ADAM.COMBSS(
  X, y, delta = nrow(X), lambda = lfamm,
  Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
  CG = CG, trunc = trunc, criterion = "FAMM", MMweights = T,
)
TimeFamm = as.numeric(Sys.time() - start, units = "secs")

start = Sys.time()
Amse = ADAM.COMBSS(
  X, y, delta = nrow(X), lambda = lmse, 
  Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
  CG = CG, trunc = trunc, criterion = "MSE"
)
TimeMSE = as.numeric(Sys.time() - start, units = "secs")

start = Sys.time()
Afrb = ADAM.COMBSS(
  X, y, delta = nrow(X), lambda = lfrb,
  Niter = Niter, alpha = alpha, epoch = epoch, tol = tol,
  CG = CG, trunc = trunc, criterion = "FRB", MMweights = T,
)
TimeFRB = as.numeric(Sys.time() - start, units = "secs")
vnames = colnames(X)
paths  = list(FAMM = Afamm, FRB = Afrb, MSE = Amse)

df = do.call(rbind, lapply(names(paths), function(nm) {
  A = paths[[nm]]
  colnames(A) = vnames
  K = nrow(A); p = ncol(A)
  prev = rbind(rep(0, p), A[-K, , drop = FALSE])
  data.frame(
    k      = rep(seq_len(K), times = p),
    var    = factor(rep(vnames, each = K), levels = vnames),
    method = nm,
    state  = as.vector(ifelse(A == 0, "Inactive",
                              ifelse(prev == 0, "Newly entered", "Active")))
  )
}))
library(ggplot2)


colnames(X)[12] = "height"
vnames = colnames(X)
paths  = list(`Robust COMBSS` = Afamm, FRB = Afrb, `COMBSS` = Amse)

df = do.call(rbind, lapply(names(paths), function(nm) {
  A = paths[[nm]]
  colnames(A) = vnames
  K = nrow(A); p = ncol(A)
  prev = rbind(rep(0, p), A[-K, , drop = FALSE])
  data.frame(
    k      = rep(seq_len(K), times = p),
    var    = factor(rep(vnames, each = K), levels = vnames),
    method = nm,
    state  = as.vector(ifelse(A == 0, "Inactive",
                              ifelse(prev == 0, "Newly Entered", "Active")))
  )
}))

df$method = factor(df$method, levels = names(paths))
df$k      = factor(df$k, levels = seq_len(max(df$k)))
df$state  = factor(df$state, levels = c("Inactive", "Active", "Newly Entered"))

pal = c("Inactive"      = "#F2F0EB",
         "Active"        = "#3C5488",
         "Newly Entered" = "#E07B54")
df$method = factor(df$method, levels = c("FRB", "Robust COMBSS", "COMBSS"))
ggplot(df, aes(var, k, fill = state)) +
  geom_tile(colour = "white", linewidth = 0.7, width = 0.95, height = 0.95) +
  facet_wrap(~ method) +
  scale_fill_manual(values = pal) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = "Subset size  k", fill = NULL) +
  coord_equal() +
  theme_few_grid(base_size = 15) +
  theme(panel.grid   = element_blank(),
        panel.border = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        strip.text   = element_text(face = "bold"),
        plot.title   = element_text(face = "bold", margin = margin(b = 10)),
        legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        panel.spacing   = unit(1.2, "lines"))

# 12x8