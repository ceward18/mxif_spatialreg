
library(ggplot2)
library(nimble)
library(grid)
library(gridExtra)
library(ggh4x)
library(knitr)
library(kableExtra)
library(dplyr)


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




# MSE, absolute difference, coverage probability, avg time
prop_stats <- readRDS('results/prop_stats.rds')

prop_stats <- merge(prop_stats, beta_vals,
                    by = 'beta_idx', all.x = T)


or_stats <- subset(prop_stats, coef == 'OR' & 
                       model_type %in% c('no_corr', 
                                         'pc_sqexp_norm_sigma',
                                         'nngp'))


or_stats$model_type <- factor(or_stats$model_type,
                              levels = c('no_corr', 
                                         'pc_sqexp_norm_sigma',
                                         'nngp'),
                              labels = c('No spatial correlation',
                                         'Eigen-decomposition',
                                         'NNGP'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))



# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p1 <- ggplot(subset(or_stats, zero_distance == 50), 
             aes(x = exp(beta_val), y = bias,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue', 'goldenrod2')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio',
         y = 'Bias',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)


# no corr vs estimation by sigma-spat, nsubject = 30, nimages = 1, zero_distance = 50
p2 <- ggplot(subset(or_stats, zero_distance == 50), 
             aes(x = exp(beta_val), y = cover,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue', 'goldenrod2')) + 
    ggtitle('Coverage of 95% CIs') +
    labs(x = 'Odds ratio',
         y = 'Coverage  probability',
         col = 'Model') +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1)


png('figures/bias_cover_single.png', units = 'in', res = 500, height =7, width = 12)
grid.arrange(p1, p2, nrow = 2)
dev.off()



time_tab <- readRDS('results/time_tab.rds')


time_tab$model_type <- factor(time_tab$model_type,
                              levels = c('no_corr', 
                                         'pc_sqexp_norm_sigma',
                                         'nngp'),
                              labels = c('No spatial correlation',
                                         'Eigen-decomposition',
                                         'NNGP'))

png('figures/comp_time.png', units = 'in', res = 500, 
    height = 3, width = 7)

ggplot(time_tab, aes(x = model_type, y = median_time,
                     ymin = lower_time, ymax = upper_time, 
                     col = model_type)) + 
    geom_point(position = position_dodge(0.9)) + 
    geom_errorbar(position = position_dodge(0.9),width = 0.3) + 
    labs(x = 'Number of subjects', y = 'Time to run (sec)',
         col = 'Model') +
    scale_color_manual(values = c('tomato', 'dodgerblue', 'goldenrod2'))  +
    theme_bw() + 
    theme(strip.text = element_text(size = 12),
          strip.background = element_rect(fill = 'white'),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 16, h = 0.5))

dev.off()
