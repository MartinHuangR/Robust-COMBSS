# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
load("COMBSS/COMBSS1C-NS.RData")
load("COMBSS/COMBSS2C-NS.RData")
source("Functions.R")

cf1 = rbind(famm1.Clean,          mse1.Clean,          frb1.Clean,
            famm1.S1,      mse1.S1,      frb1.S1,
            famm1.S2,      mse1.S2,      frb1.S2,
            famm1.S3,     mse1.S3,     frb1.S3,
            famm1.S4,     mse1.S4,     frb1.S4,
            famm1.S5, mse1.S5, frb1.S5) |> data.frame()
colnames(cf1) = paste0("Size ", 1:ncol(cf1))
cf1$Contamination = rep(c("Clean", "S1", "S2", "S3",
                          "S4", "S5"), each = 3)
cf1 = rownames_to_column(cf1, var = "tmp")
cf1$Method = "COMBSS"
cf1$Method[grepl("famm", cf1$tmp)] = "Robust COMBSS"
cf1$Method[grepl("frb", cf1$tmp)] = "FRB"
cf1$Method = factor(cf1$Method, levels = c("FRB", "Robust COMBSS", "COMBSS"))
mcf1 = reshape2::melt(cf1 |> dplyr::select(-tmp), id.vars = c("Contamination", "Method"))
mcf1$Setting = "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3"
cf2 = rbind(famm2.Clean,          mse2.Clean,          frb2.Clean,
            famm2.S1,      mse2.S1,      frb2.S1,
            famm2.S2,      mse2.S2,      frb2.S2,
            famm2.S3,     mse2.S3,     frb2.S3,
            famm2.S4,     mse2.S4,     frb2.S4,
            famm2.S5, mse2.S5, frb2.S5) |> data.frame()
colnames(cf2) = paste0("Size ", 1:ncol(cf2))
cf2$Contamination = rep(c("Clean", "S1", "S2", "S3",
                          "S4", "S5"), each = 3)
cf2 = rownames_to_column(cf2, var = "tmp")
cf2$Method = "COMBSS"
cf2$Method[grepl("famm", cf2$tmp)] = "Robust COMBSS"
cf2$Method[grepl("frb", cf2$tmp)] = "FRB"
cf2$Method = factor(cf2$Method, levels = c("FRB", "Robust COMBSS", "COMBSS"))
mcf2 = reshape2::melt(cf2 |> dplyr::select(-tmp), id.vars = c("Contamination", "Method"))

mcf2$Setting = "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5"

mcf = rbind(mcf1, mcf2)
# mcf = mcf1

mcf$size = readr::parse_number(as.character(mcf$variable)) |> as.numeric()
mcf$Contamination = factor(mcf$Contamination,
                           levels = c("Clean", "S1", "S2", "S3",
                                      "S4", "S5"),
                           labels = c("plain(Clean)",
                                      "atop(plain(S1)*':'~~c==4*','~d==1,tau==0.1)",
                                      "atop(plain(S2)*':'~~c==4*','~d==1,tau==0.2)",
                                      "atop(plain(S3)*':'~~c==4*','~d==0.1,tau==0.1)",
                                      "atop(plain(S4)*':'~~c==4*','~d==0.1,tau==0.2)",
                                      "atop(plain(S5)*':'~~c==6*','~d==0.1,tau==0.3)"))

vlines <- data.frame(
  Setting = c("(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3",
              "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5"),
  xint    = c(3, 5)
)

p1 = ggplot(mcf, aes(x = size, y = value/1000, colour = Method, group = Method)) +
  geom_line(linewidth = 0.8) +
  facet_grid(Contamination~Setting,
              labeller = labeller(Setting = label_parsed, Contamination = label_parsed), scales = "free") +
  facetted_pos_scales(
    x = list(
      Setting == "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5" ~ scale_x_continuous(breaks = c(1,5,25,49),
                                                                                      limits = c(1,49)),
      Setting == "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3" ~ scale_x_continuous(breaks = c(1,3,10,19),
                                                                                     limits = c(1,19))
    )) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ylab("Good/Correct Model Empirical Selection Probability") +
  xlab("Model Size") + 
  geom_vline(data = vlines, aes(xintercept = xint), linetype = "dashed", alpha = 0.5) + 
  theme_few_grid(base_size = 15) + 
  scale_colour_brewer(palette = "Set2")
#12x10
p1





cft1 = c(famm1.Cleant, famm1.S1t, famm1.S2t,
         famm1.S3t, famm1.S4t, famm1.S5t,
         mse1.Cleant, mse1.S1t, mse1.S2t,
         mse1.S3t, mse1.S4t, mse1.S5t,
         frb1.Cleant, frb1.S1t, frb1.S2t,
         frb1.S3t, frb1.S4t, frb1.S5t)
cft1 = data.frame(Time = cft1, method = rep(c("Robust COMBSS", "COMBSS", "FRB"), each = 6000),
                  Setting = "(i):~n==200*`,`~p==20*`,`~`|`*beta[S]*`|`==3")
cft1$method = factor(cft1$method, levels = c("FRB", "Robust COMBSS", "COMBSS"))



cft2 = c(famm2.Cleant, famm2.S1t, famm2.S2t,
         famm2.S3t, famm2.S4t, famm2.S5t,
         mse2.Cleant, mse2.S1t, mse2.S2t,
         mse2.S3t, mse2.S4t, mse2.S5t,
         frb2.Cleant, frb2.S1t, frb2.S2t,
         frb2.S3t, frb2.S4t, frb2.S5t)
cft2 = data.frame(Time = cft2, method = rep(c("Robust COMBSS", "COMBSS", "FRB"), each = 6000),
                  Setting = "(ii):~n==500*`,`~p==50*`,`~`|`*beta[S]*`|`==5")
cft2$method = factor(cft2$method, levels = c("FRB", "Robust COMBSS", "COMBSS"))

cft = rbind(cft1, cft2)
p2 = ggplot(cft, aes(x = method,
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
#5x10
p2

cft |> group_by(method) |> summarise(m = mean(Time))

