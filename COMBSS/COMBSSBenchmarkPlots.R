# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R-code functions to produce Section 4.2 Results                   #
# in the paper:                                                     #
# Robust Best Subset Selection via Fast Approximate MM-Estimation   #
#                                                                   #
# Author: Martin Huang (martin.huang@sydney.edu.au)                 #          
#         School of Mathematics & Statistics, University of Sydney  #          
#         AUSTRALIA                                                 #          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
load("COMBSS/Cb.RData")
set2_cols = RColorBrewer::brewer.pal(4, "Set2")
cols = c(
  `Robust COMBSS` = set2_cols[2],
  COMBSS  = set2_cols[3]
)
rld0$df$Tau = "Contamination: 0%"
rld0.1$df$Tau = "Contamination: 10%"
rld0.2$df$Tau = "Contamination: 20%"
rld = rbind(rld0$df,rld0.1$df,rld0.2$df)
rld$Tau = factor(rld$Tau, levels = c("Contamination: 0%","Contamination: 10%","Contamination: 20%"))
trld = mean(as.vector(rld0[[2]])/as.vector(rld0[[1]]), as.vector(rld0.1[[2]])/as.vector(rld0.1[[1]]), as.vector(rld0.2[[2]])/as.vector(rld0.2[[1]]))

ggplot(rld, aes(x = log(grid.comb))) +
  geom_line(aes(y = value, colour = method)) +
  geom_point(aes(y = value, colour = method), size = 2) +
  geom_step(aes(y = nsel), colour = "grey50", linetype = 2) +
  facet_grid(method ~ Tau, labeller = labeller(Tau = label_value)) +
  scale_x_continuous(
    name = expression(log(lambda)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  labs(y = "MSE (validation set) and Number of Selected Variables") +
  theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols) +
  theme(axis.ticks.x.top = element_blank(), legend.position = "none")


rhd0$df$Tau = "Contamination: 0%"
rhd0.1$df$Tau = "Contamination: 10%"
rhd0.2$df$Tau = "Contamination: 20%"
rhd = rbind(rhd0$df,rhd0.1$df,rhd0.2$df)
rhd$Tau = factor(rhd$Tau, levels = c("Contamination: 0%","Contamination: 10%","Contamination: 20%"))

trhd = mean(as.vector(rhd0[[2]])/as.vector(rhd0[[1]]), as.vector(rhd0.1[[2]])/as.vector(rhd0.1[[1]]), as.vector(rhd0.2[[2]])/as.vector(rhd0.2[[1]]))

ggplot(rhd, aes(x = log(grid.comb))) +
  geom_line(aes(y = value, colour = method)) +
  geom_point(aes(y = value, colour = method), size = 2) +
  geom_step(aes(y = nsel), colour = "grey50", linetype = 2) +
  facet_grid(method ~ Tau, labeller = labeller(Tau = label_value)) +
  scale_x_continuous(
    name = expression(log(lambda)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  labs(y = "MSE (validation set) and Number of Selected Variables") +
  theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols) +
  theme(axis.ticks.x.top = element_blank(), legend.position = "none")
# 16x10