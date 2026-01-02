library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(ggpubr)
library(latex2exp)

theme_set(theme_bw())
options(ggplot2.discrete.colour= c("#E69F00","#56B4E9","#CC79A7", "#009E73", "#F0E442", "#0072B2",
                                   "#D55E00","#999999", "#000000"))

########### Mean plot ################

mean_aucs <- read_csv('simulations/results/mean_multiway_aucs.csv')

plot_mean_aucs <- mean_aucs %>%
  group_by(Method, Mean) %>%
  summarize(AUC = mean(AUC)) 

mean_plot <- ggplot(plot_mean_aucs, aes(x = Mean, AUC, colour = Method)) + 
  geom_line() +
  ylim(0,1) +
  labs(x = 'Mean of C*') +
#  geom_vline(xintercept = 0.25, linetype = 'dashed', color = 'grey')
theme(legend.position = 'bottom') 

legend <- get_legend(mean_plot)
## Add vertical lines to legend



mean_plot <- mean_plot + theme(legend.position = 'none')

########### Rho time plot ################
rho_aucs <- read_csv('simulations/results/rho_multiway_aucs.csv')

plot_rho_aucs  <- rho_aucs %>%
  group_by(Method, Rho) %>%
  summarize(AUC = mean(AUC)) 

rho_plot <- ggplot(plot_rho_aucs, aes(x = Rho, AUC, colour = Method)) + 
  geom_line()+
  labs(x = expression(rho[t] ~ ': corr of ' ~ beta^(1) ~ ','~ beta^(2) )) +
  theme(legend.position = 'none') +
  ylim(0,1) +
  geom_vline(xintercept = 0.8, linetype = 'dashed', color = 'grey') + 
  ## smooth lines
  geom_smooth(se=FALSE, method = 'lm')


########### Cluster plot ################
cluster_aucs <- read_csv('simulations/results/cluster_multiway_aucs.csv')

plot_cluster_aucs  <- cluster_aucs %>%
  group_by(Method, Cluster) %>%
  summarize(AUC = mean(AUC)) 

cluster_plot <- ggplot(plot_cluster_aucs, aes(x = Cluster, AUC, colour = Method)) + 
  geom_line()+
  labs(x = 'Number of Clusters') +
  theme(legend.position = 'none') +
  ylim(0,1) + 
  geom_vline(xintercept = 7, linetype = 'dashed', color = 'grey') + 
  geom_smooth(se=FALSE, method = 'lm')


combined <- grid.arrange(mean_plot, rho_plot, cluster_plot, legend, nrow = 1,widths=c(3,3,3,1.2))


g <- combined#grid.arrange(combined, , nrow = 1, heights = c(10,1.5))
ggsave('figures/supplementary_1-test.png',g, width = 8, height = 5, units = 'in')
