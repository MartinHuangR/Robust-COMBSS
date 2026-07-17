# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
load("COMBSS/COMBSS1.RData")
load("COMBSS/COMBSS2.RData")
source("Functions.R")

cf1 = rbind(famm1.0,    mse1.0,    frb1.0,
            famm1.0.1,  mse1.0.1,  frb1.0.1) |> data.frame()
colnames(cf1) = paste0("Size ", 1:ncol(cf1))
cf1$Tau = rep(c(0,0.1), each = 3)
cf1$klev = 2
cf1$kslo = 100
cf1 = rownames_to_column(cf1, var = "tmp")
cf1$Method = "MSE"
cf1$Method[grepl("famm", cf1$tmp)] = "FAMM"
cf1$Method[grepl("frb", cf1$tmp)] = "FRB"
cf1$Method = factor(cf1$Method, levels = c("FRB", "FAMM", "MSE"))
mcf1 = reshape2::melt(cf1 |> dplyr::select(-tmp), id.vars = c("Tau",
                                                              "kslo", "klev", "Method"))
mcf1$Setting = "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3"
cf2 = rbind(famm2.0,    mse2.0,    frb2.0,
            famm2.0.1,  mse2.0.1,  frb2.0.1) |> data.frame()
colnames(cf2) = paste0("Size ", 1:ncol(cf2))
cf2$Tau = rep(c(0,0.1), each = 3)
cf2$klev = 2
cf2$kslo = 100
cf2 = rownames_to_column(cf2, var = "tmp")
cf2$Method = "MSE"
cf2$Method[grepl("famm", cf2$tmp)] = "FAMM"
cf2$Method[grepl("frb", cf2$tmp)] = "FRB"
cf2$Method = factor(cf2$Method, levels = c("FRB", "FAMM", "MSE"))
mcf2 = reshape2::melt(cf2 |> dplyr::select(-tmp), id.vars = c("Tau",
                                                              "kslo", "klev", "Method"))

mcf2$Setting = "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5"

mcf = rbind(mcf1, mcf2)


mcf$size = readr::parse_number(as.character(mcf$variable)) |> as.numeric()
mcf$TauP = paste0("Contamination: ", mcf$Tau * 100, "%")
mcf$TauP = factor(mcf$TauP, levels = c("Contamination: 0%", "Contamination: 10%"))
  
vlines <- data.frame(
  Setting = c("(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3",
              "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5"),
  xint    = c(3, 5)
)
  
ggplot(mcf, aes(x = size, y = value/1000, colour = Method, group = Method)) +
  geom_line(linewidth = 0.8) +
  facet_grid(TauP~Setting,labeller = labeller(Setting = label_parsed, TauP = label_value), scales = "free") +
  facetted_pos_scales(
    y = list(
      TauP == "Contamination: 0%"   ~ scale_y_continuous(limits = c(0.5, 1), breaks = c(0.5,0.6,0.7,0.8,0.9,1)),
      TauP == "Contamination: 10%" ~ scale_y_continuous(limits = c(0.5, 1))),
    x = list(
      Setting == "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5" ~ scale_x_continuous(breaks = c(1,5,25,49),
                                    limits = c(1,49)),
      Setting == "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3" ~ scale_x_continuous(breaks = c(1,3,10,19),
                                                                                 limits = c(1,19))
    )) + 
  ylab("Correct Model Empirical Selection Probability") +
  xlab("Model Size") + 
  geom_vline(data = vlines, aes(xintercept = xint), linetype = "dashed", alpha = 0.5) + 
  theme_few_grid(base_size = 15) + 
  scale_colour_brewer(palette = "Set2")





cft1 = c(famm1.0t, famm1.0.1t,  
         mse1.0t,  mse1.0.1t, 
         frb1.0t,  frb1.0.1t)
cft1 = data.frame(Time = cft1, method = rep(c("FAMM", "MSE", "FRB"), each = 2000),
                  Setting = "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3")
cft1$method = factor(cft1$method, levels = c("FRB", "FAMM", "MSE"))



cft2 = c(famm2.0t, famm2.0.1t,  
         mse2.0t,  mse2.0.1t, 
         frb2.0t,  frb2.0.1t)
cft2 = data.frame(Time = cft2, method = rep(c("FAMM", "MSE", "FRB"), each = 2000),
                  Setting = "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5")
cft2$method = factor(cft2$method, levels = c("FRB", "FAMM", "MSE"))

cft = rbind(cft1, cft2)
ggplot(cft, aes(x = method,
                       y = log10(as.numeric(Time)),
                       colour = method)) +
  geom_boxplot() +
  facet_grid(~Setting, labeller = labeller(Setting = label_parsed)) + 
  ylab("Seconds to run 1 dataset (Log scale)") +
  xlab("Method") +
  theme_few_grid(base_size = 15) +
  scale_colour_brewer(palette = "Set2") + 
  scale_y_continuous(
    breaks = log10(c(1, 10, 100, 1000)),
    labels = c("1", "10", "100", "1000")
  ) +
  theme(legend.position = "none")
