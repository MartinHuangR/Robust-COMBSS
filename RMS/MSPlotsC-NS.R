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
load("RMS/MSC-NS.RData")
set2_cols = RColorBrewer::brewer.pal(4, "Set2")
cols = c(
  FRB  = set2_cols[1],
  FAMM = set2_cols[2],
  MM  = "#7C5A3C"
)


SC1 = rbind(
  FAMMt0Cleann500UO, FRBt0Cleann500UO,
  FAMMt0.1S1n500UO, FRBt0.1S1n500UO,
  FAMMt0.2S2n500UO, FRBt0.2S2n500UO,
  FAMMt0.1S3n500UO, FRBt0.1S3n500UO,
  FAMMt0.2S4n500UO, FRBt0.2S4n500UO,
  FAMMt0.3S5n500UO, FRBt0.3S5n500UO,
  MMt0Cleann500UO,
  MMt0.1S1n500UO,
  MMt0.2S2n500UO,
  MMt0.1S3n500UO,
  MMt0.2S4n500UO,
  MMt0.3S5n500UO)


SC2 = rbind( FAMMt0Cleann1000UO, FRBt0Cleann1000UO,
             FAMMt0.1S1n1000UO, FRBt0.1S1n1000UO,
             FAMMt0.2S2n1000UO, FRBt0.2S2n1000UO,
             FAMMt0.1S3n1000UO, FRBt0.1S3n1000UO,
             FAMMt0.2S4n1000UO, FRBt0.2S4n1000UO,
             FAMMt0.3S5n1000UO, FRBt0.3S5n1000UO,
             MMt0Cleann1000UO,
             MMt0.1S1n1000UO,
             MMt0.2S2n1000UO,
             MMt0.1S3n1000UO,
             MMt0.2S4n1000UO,
             MMt0.3S5n1000UO)

# Plotmath labels: a is shift, b is spread, and tau is the contamination proportion.
SC1$Scenario = ifelse(SC1$Tau == 0, "plain(Clean)",
                      paste0("c==", SC1$shift, "*','~d==", SC1$spread,
                             "*','~tau==", SC1$Tau))
SC1$Scenario = factor(SC1$Scenario, levels = unique(SC1$Scenario))
SC2$Scenario = ifelse(SC2$Tau == 0, "plain(Clean)",
                      paste0("c==", SC2$shift, "*','~d==", SC2$spread,
                             "*','~tau==", SC2$Tau))
SC2$Scenario = factor(SC2$Scenario, levels = unique(SC2$Scenario))

# Short horizontal ticks; give the full settings beside each legend shape.
scenario_ticks = c("Clean", "S1", "S2", "S3", "S4", "S5")
names(scenario_ticks) = levels(SC1$Scenario)
scenario_legend = c("plain(Clean)",
                    paste0("plain(", scenario_ticks[-1], ")*':'~~",
                           levels(SC1$Scenario)[-1]))
scenario_legend = parse(text = scenario_legend)
names(scenario_legend) = levels(SC1$Scenario)

SC1$Method = factor(SC1$method, levels = c("FRB", "FAMM", "MM"))
SC1$Tau = factor(SC1$Tau)
SC1$Criterion <- factor(
  SC1$Criterion,
  levels = c("M02.ind", "M05.ind"),
  labels = c("~tilde(alpha)[250*','*500]", "~hat(alpha)[250*','*500]")
)
SC1 = SC1[SC1$Criterion == "~hat(alpha)[250*','*500]",]
SC1p1 = ggplot(SC1, aes(x = Scenario, y = True/1000, colour = Method, group = Method)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(labels = scenario_ticks) +
  ylim(0,1) + 
  theme_bw() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Contamination Setting") +
  theme_few_grid(base_size = 15) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))+ 
  scale_colour_manual(values = cols)
SC1p2 = ggplot(SC1, aes(x = as.numeric(TotalTime)/60, y = True/1000, colour = Method, shape = Scenario)) + geom_point(size = 3) +
  scale_shape_discrete(labels = scenario_legend) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, ncol = 2, byrow = T)) +
  ylim(0,1) + 
  scale_x_log10(breaks = c(5,10,20,100,150), limits = c(5,150)) +
  ylab("True Model Empirical Selection Probability") + xlab("Minutes taken to run 1000 datasets") +
  theme(legend.position = "bottom") + theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols)

SC2$Method = factor(SC2$method, levels = c("FRB", "FAMM","MM"))
SC2$Tau = factor(SC2$Tau)
SC2$Criterion <- factor(
  SC2$Criterion,
  levels = c("M02.ind", "M05.ind"),
  labels = c("~tilde(alpha)[250*','*500]", "~hat(alpha)[250*','*500]")
)
SC2 = SC2[SC2$Criterion == "~hat(alpha)[250*','*500]",]
SC2p1 = ggplot(SC2, aes(x = Scenario, y = True/1000, colour = Method, group = Method)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(labels = scenario_ticks) +
  theme_bw() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Contamination Setting") +
  theme_few_grid(base_size = 15) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))+ 
  scale_colour_manual(values = cols)
SC2p2 = ggplot(SC2, aes(x = as.numeric(TotalTime)/60, y = True/1000, colour = Method, shape = Scenario)) + geom_point(size = 3) +
  scale_shape_discrete(labels = scenario_legend) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, ncol = 2, byrow = T)) +
  ylim(0,1) + 
  scale_x_log10(breaks = c(20,30,50,300,500), limits = c(20,500)) +
  ylab("True Model Empirical Selection Probability") + xlab("Minutes taken to run 1000 datasets") +
  theme(legend.position = "bottom") + theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols)

design <- "
AB
CD
GG
"

SC1p1_noleg <- SC1p1 +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none", shape = "none", linetype = "none")

SC2p1_noleg <- SC2p1 +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none", shape = "none", linetype = "none")

p <-
  (
    SC1p1_noleg + SC2p1_noleg +
      SC1p2 + SC2p2 +
      guide_area()
  ) +
  plot_layout(
    design = design,
    guides = "collect",
    axis_titles = "collect",
    widths = c(1, 1),
    heights = c(1, 1, 0.6)
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.text = element_text(size = 15),
    plot.tag.position = "topleft",
    plot.tag = element_text(hjust = 0, vjust = 1, size = 20)
  )

# 8 x 10
fig2 = p + plot_annotation(tag_levels = list(c("(i)", "(ii)", "", "")))
fig2
# ggsave("FAMM_paper/Figures/SC4.pdf", fig2, width = 8, height = 10)
