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

set2_cols = RColorBrewer::brewer.pal(4, "Set2")
cols = c(
  FRB  = set2_cols[1],
  FAMM = set2_cols[2],
  MM  = "#7C5A3C"
)


load("RMS/MSsm.RData")
SM = data.frame(rbind(FRB, FAMM,MM),
                Method = rep(c("FRB", "FAMM", "MM"), each = nrow(FRB)))
SM$Method = rep(c("FRB", "FAMM", "MM"), each = nrow(FRB))
SM$Method = factor(SM$Method, levels = c("FRB", "FAMM", "MM"))
SM$Outlier = factor(SM$Error, levels = c("N", "1/4", "1/8", "3/8", "C", "S"))
SM$Criterion <- factor(
  SM$Criterion,
  levels = c("M02.ind", "M05.ind"),
  labels = c("~tilde(alpha)[24*','*64]", "~hat(alpha)[24*','*64]")
)
SM = SM[SM$Criterion == "~hat(alpha)[24*','*64]",]
SM$Error = factor(SM$Error, levels = c("N","1/8", "1/4", "3/8",  "C", "S"))
p1 = ggplot(SM, aes(x = Error, y = True/1000, colour = Method, group = Method)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Outlier Function") +
  theme_few_grid(base_size = 15) +
  theme(legend.position = "none")+ 
  ylim(0.2,1) + 
  scale_colour_manual(values = cols)
p2 = ggplot(SM, aes(x = as.numeric(TotalTime), y = True/1000, colour = Method, shape = Error)) +
  geom_point(size = 3) +
  theme_bw() +
  ylim(0.2, 1) +
  scale_x_log10() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Seconds taken to run 1000 datasets") +
  theme(legend.position = "bottom") + theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols)

# 5x11
fig1 = p1 + p2 + plot_layout(guides = 'collect', axis_title = "collect")
fig1
# ggsave("FAMM_paper/Figures/SM1.pdf", fig1, width = 12, height = 5)

load("RMS/MS1.RData")
S1 = rbind(FAMMt0s1UO, FRBt0s1UO, MMt0s1UO,
           FAMMt0.05s1UO, FRBt0.05s1UO, MMt0.05s1UO,
           FAMMt0.1s1UO, FRBt0.1s1UO, MMt0.1s1UO,
           FAMMt0.15s1UO, FRBt0.15s1UO, MMt0.15s1UO,
           FAMMt0.2s1UO, FRBt0.2s1UO, MMt0.2s1UO)
S1$Method = factor(S1$method, levels = c("FRB", "FAMM", "MM"))
S1$Tau = factor(S1$Tau)
S1$Criterion <- factor(
  S1$Criterion,
  levels = c("M02.ind", "M05.ind"),
  labels = c("~tilde(alpha)[250*','*500]", "~hat(alpha)[250*','*500]")
)
S1 = S1[S1$Criterion == "~hat(alpha)[250*','*500]",]
s1p1 = ggplot(S1, aes(x = Tau, y = True/1000, colour = Method, group = Method)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Outlier Proportion") +
  theme_few_grid(base_size = 15) +
  theme(legend.position = "none")+ 
  ylim(0.97,1) + 
  scale_colour_manual(values = cols)
s1p2 = ggplot(S1, aes(x = as.numeric(TotalTime)/60, y = True/1000, colour = Method, shape = Tau)) + geom_point(size = 3) +
  ylim(0.97,1) + 
  scale_x_log10(breaks = c(10,30,100,200),limits = c(8, 200)) + 
  ylab("True Model Empirical Selection Probability") + xlab("Minutes taken to run 1000 datasets") +
  theme(legend.position = "bottom") + theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols)


load("RMS/MS2.RData")
S2 = rbind(FAMMt0s2UO, FRBt0s2UO, MMt0s2UO,
           FAMMt0.05s2UO, FRBt0.05s2UO, MMt0.05s2UO,
           FAMMt0.1s2UO, FRBt0.1s2UO, MMt0.1s2UO,
           FAMMt0.15s2UO, FRBt0.15s2UO, MMt0.15s2UO,
           FAMMt0.2s2UO, FRBt0.2s2UO, MMt0.2s2UO)

S2$Method = factor(S2$method, levels = c("FRB", "FAMM", "MM"))
S2$Tau = factor(S2$Tau)
S2$Criterion <- factor(
  S2$Criterion,
  levels = c("M02.ind", "M05.ind"),
  labels = c("~tilde(alpha)[500*','*1000]", "~hat(alpha)[500*','*1000]")
)
S2 = S2[S2$Criterion == "~hat(alpha)[500*','*1000]",]
s2p1 = ggplot(S2, aes(x = Tau, y = True/1000, colour = Method, group = Method)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("True Model Empirical Selection Probability") +
  xlab("Outlier Proportion") +
  theme_few_grid(base_size = 15) +
  theme(legend.position = "none")+ 
  ylim(0.97,1) + 
  scale_colour_manual(values = cols)

s2p2 = ggplot(S2, aes(x = as.numeric(TotalTime)/60, y = True/1000, colour = Method, shape = Tau)) + geom_point(size = 3) +
  theme_bw() +
  scale_x_log10(breaks = c(30, 100, 500,1000), limits = c(30, 1000)) +
  ylim(0.97,1) + 
  ylab("True Model Empirical Selection Probability") + xlab("Minutes taken to run 1000 datasets") +
  theme(legend.position = "bottom") + theme_few_grid(base_size = 15) +
  scale_colour_manual(values = cols)

design <- "
ABG
CDG
"

s1p1_noleg <- s1p1 +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none", shape = "none", linetype = "none")

s2p1_noleg <- s2p1 +
  theme(legend.position = "none") +
  guides(fill = "none", color = "none", shape = "none", linetype = "none")

p <-
  (
    s1p1_noleg + s2p1_noleg +
      s1p2 + s2p2 +
      guide_area()
  ) +
  plot_layout(
    design = design,
    guides = "collect",
    axis_titles = "collect",
    widths = c(1, 1, 0.28)
  ) &
  theme(
    legend.position = "right",
    plot.tag.position = "topleft",
    plot.tag = element_text(hjust = 0, vjust = 1, size = 20)
  )

# 8x10
fig2 = p + plot_annotation(tag_levels = list(c("(i)", "(ii)", "", "")))
fig2
# ggsave("FAMM_paper/Figures/S14.pdf", fig2, width = 10, height = 8)
