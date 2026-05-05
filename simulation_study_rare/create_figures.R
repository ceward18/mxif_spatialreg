################################################################################
# create figures for presentation
################################################################################

library(ggplot2)
library(nimble)
library(grid)
library(gridExtra)
library(ggh4x)
library(knitr)
library(kableExtra)
library(dplyr)
library(tidyverse)
library(parallel)
library(foreach)
library(MASS)
library(spatstat.geom)
library(spatstat.random)
library(mvnfast)
library(Matrix)
library(doSNOW)
library(cowplot)

source('sim_function.R')
source('helper_functions.R')

################################################################################
# Figure S4 - bias and coverage of 95% CI

# prop outcome, beta1 corresponds to OR
ors <- c(1, 1.25, 1.5, 2, 3)
ors <- c(rev(1/ors[-1]), ors)

beta_vals <- data.frame(beta_idx = 1:length(ors),
                        beta_val = log(ors))


# MSE, absolute difference, coverage probability, avg time
prop_stats <- readRDS('results/prop_stats.rds')


prop_stats$n_subjects <- factor(prop_stats$n_subjects)
prop_stats <- merge(prop_stats, beta_vals,
                    by = 'beta_idx', all.x = T)


or_stats <- subset(prop_stats, 
                   model_type %in% c('no_corr', 
                                     'pc_sqexp',
                                     'inla'))


or_stats$model_type <- factor(or_stats$model_type,
                              levels = c('inla',
                                         'pc_sqexp', 
                                         'no_corr'),
                              labels = c('INLA-SPDE',
                                         'Eigen-decomposition',
                                         'No spatial correlation'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))


fig2_theme <- theme_bw() + 
    theme(strip.text = element_text(size = 16),
          strip.background = element_rect(fill = 'white'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 17),
          legend.key.width = unit(1, "cm"),   
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.title = element_text(size = 18), 
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 20, h = 0.5))


# relative Bias 
p1 <- ggplot(subset(or_stats, coef == 'OR_R'), 
             aes(x = exp(beta_val), y = bias_rel,  col = model_type)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1) + 
    # geom_point(size = 3, alpha = 0.8) + 
    geom_line(linewidth = 1.5, alpha = 0.8) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'dodgerblue', 'tomato')) + 
    ggtitle('Relative Bias') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Relative Bias',
         col = 'Model') +
    fig2_theme + 
    scale_y_continuous(limits = c(-0.4, 0.4))


# Bias 
p2 <- ggplot(subset(or_stats, coef == 'OR_R'), 
             aes(x = exp(beta_val), y = bias,  col = model_type)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1) + 
    # geom_point(size = 3, alpha = 0.8) + 
    geom_line(linewidth = 1.5, alpha = 0.8) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'dodgerblue', 'tomato')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Bias',
         col = 'Model') +
    fig2_theme + 
    scale_y_continuous(limits = c(-0.8, 0.8))

# coverage
p3 <- ggplot(subset(or_stats, coef == 'OR_R' ), 
             aes(x = exp(beta_val), y = cover,  col = model_type)) +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1) +
    # geom_point(size = 3, alpha = 0.8) + 
    geom_line(linewidth = 1.5, alpha = 0.8) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'dodgerblue', 'tomato')) + 
    ggtitle('Coverage of 95% CIs') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Coverage  probability',
         col = 'Model')+ 
    fig2_theme + 
    scale_y_continuous(
        breaks = sort(unique(c(seq(0, 1, by = 0.2), 0.95))),
        labels = function(x) ifelse(x == 0.95, "0.95", x)
    )

# Extract the legend from one plot
legend <- get_legend(
    p1 + theme(legend.position = "right")
)


# Remove legends from the individual panels
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")
p3_noleg <- p3 + theme(legend.position = "none")

# Combine plots vertically, add shared legend to the right
fig_s4 <- plot_grid(
    plot_grid(p1_noleg, p2_noleg, p3_noleg, ncol = 1, align = "v"),
    legend,
    rel_widths = c(1, 0.3)
)


ggsave("figures/fig_s4_rare.png", plot = fig_s4, device = "png", 
       width = 12, height = 12, units = "in")


