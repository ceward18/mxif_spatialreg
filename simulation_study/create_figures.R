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


res_all <- readRDS('results/res_all.rds')

col_pal <- c('red', 'dodgerblue', 'royalblue', 'lightgoldenrod', 'gold4')


theme_set(theme_bw() + 
              theme(strip.text = element_text(size = 15),
                    strip.background = element_rect(fill = 'white'),
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 15),
                    plot.title = element_text(size = 16, h = 0.5)))

# prop outcome, beta1 corresponds to OR
ors <- c(1, 1.1, 1.25, 1.5, 2, 3)
ors <- c(rev(1/ors[-1]), ors)

beta_vals <- data.frame(beta_idx = 1:11,
                        beta_val = log(ors))



sigma_est_tab <- subset(res_prop, coef == 'OR') %>%
    group_by(sigma_spat, model_type) %>%
    summarise(nsims = length(sim_number),
              mse = mean((sigma_spat_est - sigma_spat)^2),
              bias = mean(sigma_spat_est - sigma_spat),
              abs_diff = mean(abs(sigma_spat_est - sigma_spat)),
              avg_est = mean(sigma_spat_est),
              truth = mean(sigma_spat)) %>%
    data.frame()
sigma_est_tab

# average number of cells
subset(res_prop, coef == 'OR' & model_type == 'no_corr') %>%
    group_by(n_subjects, n_image_sub) %>%
    summarise(avg_n_cells = median(n_cells)) %>%
    data.frame()


# MSE, absolute difference, coverage probability, avg time
prop_stats <- res_all %>%
    group_by(n_subjects, n_image_sub, zero_distance, sigma_spat, coef, model_type, beta_idx) %>%
    summarise(nsims = length(sim_number),
              mse = mean((est/truth - truth/truth)^2),
              bias = mean((est - truth)/truth),
              abs_diff = mean(abs(est/truth - truth/truth)),
              cover = mean(truth > lower & truth < upper),
              time = mean(time),
              truth = mean(truth)) %>%
    data.frame()

prop_stats$n_subjects <- factor(prop_stats$n_subjects)
prop_stats <- merge(prop_stats, beta_vals,
                    by = 'beta_idx', all.x = T)

ggplot(subset(prop_stats, coef == 'OR' & n_subjects == 10 & 
                  model_type %in% c('no_corr', 
                                    'pc_sqexp_norm_overall_sigma',
                                    'pc_sqexp_norm_overall_1',
                                    'pc_sqexp_norm_image_sigma') & 
                  n_image_sub == 1 & zero_distance == 50), 
       aes(x = exp(beta_val), y = mse,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ factor(sigma_spat ) ) + 
    ggtitle('Absolute Bias') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Bias relative to true OR',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)

ggplot(subset(prop_stats, coef == 'OR' & n_subjects == 10 & 
                  n_image_sub == 1 & zero_distance == 50), 
       aes(x = exp(beta_val), y = mse,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ factor(sigma_spat ) ) + 
    ggtitle('Absolute Bias') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Bias relative to true OR',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)


or_stats <- subset(prop_stats, coef == 'OR' & model_type %in% c('no_corr', 
                                                                'pc_est_sqexp_sigma_weight'))


or_stats$model_type <- factor(or_stats$model_type,
                              labels = c('No spatial correlation',
                                         'Eigen-decomposition offset'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))



# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p1 <- ggplot(subset(or_stats, n_subjects == 50 & 
                        n_image_sub == 5 & zero_distance == 50), 
             aes(x = exp(beta_val), y = abs(bias),  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue')) + 
    ggtitle('Absolute Bias') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Bias relative to true OR',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)


# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p2 <- ggplot(subset(or_stats, n_subjects == 50 & 
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

png('figures/mse_cover.png', units = 'in', res = 500, height =7, width = 12)
grid.arrange(p1, p2, nrow = 2)
dev.off()

or_stats$n_image_sub <- factor(or_stats$n_image_sub, 
                               labels = c('1 image per subject',
                                          '5 images per subject'))


png('figures/mse_datasize.png', units = 'in', res = 500, height =9, width = 9)
ggplot(subset(or_stats, zero_distance == 50 &
                  model_type %in% c('Eigen-decomposition offset')), 
       aes(x = exp(beta_val), y = mse,  linetype = n_subjects)) +
    # geom_point(size = 2) +
    geom_line(linewidth = 1) +
    facet_nested(sigma_spat_fac ~ n_image_sub ) +
    scale_color_manual(values = col_pal) + 
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'MSE relative to true OR',
         linetype = 'Number of subjects') +
    geom_hline(yintercept = 0) +
    theme(legend.key.width = unit(2,"cm"))
dev.off()


sigma_est_tab <- subset(res_prop, coef == 'OR' & model_type %in% c('no_corr', 'pc_est_sqexp')) %>%
    group_by(model_type, n_subjects, n_image_sub) %>%
    summarise(avg_time = mean(time),
              lower_time = quantile(time, 0.25),
              upper_time = quantile(time, 0.75)) %>%
    data.frame()
sigma_est_tab

sigma_est_tab$n_image_sub <- factor(sigma_est_tab$n_image_sub, 
                                    labels = c('1 image per subject',
                                               '5 images per subject'))


sigma_est_tab$model_type <- factor(sigma_est_tab$model_type,
                                   labels = c('No spatial correlation',
                                              'Eigen-decomposition offset'))

png('figures/comp_time.png', units = 'in', res = 500, 
    height = 3.5, width = 7)

ggplot(sigma_est_tab, aes(x = n_subjects, y = avg_time,
                          ymin = lower_time, ymax = upper_time, 
                          col = model_type)) + 
    geom_point(position = position_dodge(0.9)) + 
    geom_errorbar(position = position_dodge(0.9)) + 
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
