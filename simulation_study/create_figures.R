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

source('sim_function.R')
source('helper_functions.R')



theme_set(theme_bw() + 
              theme(strip.text = element_text(size = 15),
                    strip.background = element_rect(fill = 'white'),
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 15),
                    plot.title = element_text(size = 16, h = 0.5)))



### example datasets

rho_100 <- optim(0, function(x) (exp(-100^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par

get_sims <- function(seed, rho) {
    set.seed(seed)
    sim_data1 <- sim_data_fn(betas = c( logit(0.05), 1.2), 
                             nSub = 1, 
                             nImagePerSub = 1, 
                             rho = rho,
                             sigma_spat = 0.5,    # spatial SD
                             sigma_sub = 0.4,            # between subjects SD
                             sigma_image = 0.4)          # between images SD
    sim_data1$spat_cor <- 'Low spatial correlation'
    set.seed(seed)
    sim_data2 <- sim_data_fn(betas = c( logit(0.05), 1.2), 
                             nSub = 1, 
                             nImagePerSub = 1, 
                             rho = rho,
                             sigma_spat = 1,    # spatial SD
                             sigma_sub = 0.4,            # between subjects SD
                             sigma_image = 0.4)          # between images SD
    sim_data2$spat_cor <- 'Medium spatial correlation'
    set.seed(seed)
    sim_data3 <- sim_data_fn(betas = c( logit(0.05), 1.2), 
                             nSub = 1, 
                             nImagePerSub = 1, 
                             rho = rho,
                             sigma_spat = 2,    # spatial SD
                             sigma_sub = 0.4,            # between subjects SD
                             sigma_image = 0.4)          # between images SD
    sim_data3$spat_cor <- 'High spatial correlation'
    
    
    sim_data <- rbind.data.frame(sim_data1, sim_data2, sim_data3)
    sim_data <- subset(sim_data, numNeighbors > 0)
    sim_data$image_spec <- seed
    sim_data$n_cells <- nrow(sim_data1)
    sim_data
}

sim_dat1 <- get_sims(1, rho_100)
sim_dat2 <- get_sims(2, rho_100)
sim_dat3 <- get_sims(3, rho_100)


sim_data_all <- rbind.data.frame(sim_dat1, sim_dat2, sim_dat3)


sim_data_all$spat_cor <- factor(sim_data_all$spat_cor, 
                            levels = c('Low spatial correlation',
                                       'Medium spatial correlation',
                                       'High spatial correlation'))
sim_data_all$image_lab <- paste0('Simulation ', sim_data_all$image_spec, '\nn cells = ', sim_data_all$n_cells)

sim_data_all <- sim_data_all[order(sim_data_all$outcome),]


png('figures/sim_datasets.png', units = 'in', res = 500, height =7, width = 10)
ggplot(sim_data_all, aes(x = x, y = y, col = outcome/numNeighbors)) + 
    geom_point(size = 1) + 
    scale_color_gradient(low = 'grey80', high = 'red') + 
    facet_grid(image_lab~spat_cor) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.text = element_text(size = 12)) +
    labs(color = "Outcome") + 
    scale_x_continuous(breaks = seq(0, 750, 250))+ 
    scale_y_continuous(breaks = seq(0, 750, 250))
dev.off()

################################################################################


# prop outcome, beta1 corresponds to OR
ors <- c(1, 1.1, 1.25, 1.5, 2, 3)
ors <- c(rev(1/ors[-1]), ors)

beta_vals <- data.frame(beta_idx = 1:11,
                        beta_val = log(ors))


# MSE, absolute difference, coverage probability, avg time
prop_stats <- readRDS('results/prop_stats.rds')


prop_stats$n_subjects <- factor(prop_stats$n_subjects)
prop_stats <- merge(prop_stats, beta_vals,
                    by = 'beta_idx', all.x = T)


or_stats <- subset(prop_stats, coef == 'OR' & 
                       model_type %in% c('no_corr', 
                                         'pc_sqexp_norm_overall_sigma'))


or_stats$model_type <- factor(or_stats$model_type,
                              labels = c('No spatial correlation',
                                         'Eigen-decomposition'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))

################################################################################
# Figure 1 - bias and coverage of 95% CIs


# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p1 <- ggplot(subset(or_stats, n_subjects == 30 & 
                        n_image_sub == 5 & zero_distance == 50), 
             aes(x = exp(beta_val), y = bias,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Bias',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)


# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p2 <- ggplot(subset(or_stats, n_subjects == 30 & 
                        n_image_sub == 5 & zero_distance == 50), 
             aes(x = exp(beta_val), y = cover,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue')) + 
    ggtitle('Coverage of 95% CIs') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Coverage  probability',
         col = 'Model') +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1)

png('figures/bias_cover.png', units = 'in', res = 500, height =7, width = 12)
grid.arrange(p1, p2, nrow = 2)
dev.off()

################################################################################
# Figure 2 - MSE by dataset size



or_stats$n_image_sub <- factor(or_stats$n_image_sub, 
                               labels = c('1 image per subject',
                                          '5 images per subject'))


png('figures/mse_datasize.png', units = 'in', res = 500, height =4, width = 12)
ggplot(subset(or_stats, zero_distance == 50 &
                  model_type %in% c('Eigen-decomposition')), 
       aes(x = exp(beta_val), y = mse_rel, col = n_image_sub,  linetype = n_subjects)) +
    # geom_point(size = 2) +
    geom_line(linewidth = 1) +
    facet_nested(  ~ sigma_spat_fac) +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'MSE relative to true OR',
         linetype = 'Number of subjects',
         col = '') +
    ylim(0, 0.2) +
    theme(legend.key.width = unit(2,"cm")) + 
    scale_color_manual(values = c(rgb(255, 204, 51, max = 255),
                                  rgb(122, 0, 25, max = 255)))
dev.off()

# time to compute

time_tab <- readRDS('results/time_tab.rds')

time_tab$median_IQR <- paste0(sprintf("%.2f", round(time_tab$median_time, 2)), 
                              ' (',
                              sprintf("%.2f", round(time_tab$lower_time, 2)),
                              ', ',
                              sprintf("%.2f", round(time_tab$upper_time, 2)),
                              ')')


time_tab <- time_tab[,c('n_subjects', 'n_image_sub', 'model_type', 'median_IQR')]

time_tab_wide <- pivot_wider(time_tab,
                             id_cols = c(n_subjects, n_image_sub),
                             names_from=model_type,
                             values_from=median_IQR)


time_tab$n_image_sub <- factor(time_tab$n_image_sub, 
                                    labels = c('1 image per subject',
                                               '5 images per subject'))


time_tab$model_type <- factor(time_tab$model_type,
                                   labels = c('No spatial correlation',
                                              'Eigen-decomposition'))

png('figures/comp_time.png', units = 'in', res = 500, 
    height = 3, width = 7)

ggplot(time_tab, aes(x = factor(n_subjects), y = median_time,
                          ymin = lower_time, ymax = upper_time, 
                          col = model_type)) + 
    geom_point(position = position_dodge(0.9)) + 
    geom_errorbar(position = position_dodge(0.9),width = 0.3) + 
    facet_wrap(~n_image_sub) +
    labs(x = 'Number of subjects', y = 'Time to run (sec)',
         col = 'Model') +
    scale_color_manual(values = c('tomato', 'dodgerblue'))  +
    theme_bw() + 
    theme(strip.text = element_text(size = 12),
          strip.background = element_rect(fill = 'white'),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 16, h = 0.5))

dev.off()
