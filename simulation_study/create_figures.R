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



theme_set(theme_bw() + 
              theme(strip.text = element_text(size = 15),
                    strip.background = element_rect(fill = 'white'),
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 15),
                    plot.title = element_text(size = 16, h = 0.5)))


################################################################################
# Figure 0 - example datasets and contribution of spatial correlation from normalized offset


rho_100 <- optim(0, function(x) (exp(-100^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par

set.seed(100)
sim_data <- sim_data_fn(betas = c( logit(0.05), 1.2), 
                        nSub = 4, 
                        nImagePerSub = 1, 
                        rho = rho_100,
                        sigma_spat = 2,    # spatial SD
                        sigma_sub = 0.4,            # between subjects SD
                        sigma_image = 0.4)          # between images SD

# do eigen-decomposition

unique_images <- unique(sim_data$imageID)

# estimate variograms for each image in parallel
my_cluster <- makeCluster(4)
doParallel::registerDoParallel(cl = my_cluster)

# for each image calculate sill and range of variogram to use to get 
# estimates of zero_distance and sigma_spat
# only needs to be estimated 

vario_ests <- foreach(j = 1:length(unique_images), 
                      .combine = 'rbind') %dopar% {
                          
                          library(sp)
                          library(gstat)
                          
                          source('helper_functions.R')
                          
                          variogram_est(fov_data = sim_data[sim_data$imageID == unique_images[j],])
                          
                      }

sigma_spat_use <- mean(vario_ests[,1])
zero_dist_use <- mean(vario_ests[,2])

spatial_offsets <- foreach(j = 1:length(unique_images), 
                           .combine = 'c') %dopar% {
                               
                               source('helper_functions.R')
                               
                               image_dat <- sim_data[sim_data$imageID == unique_images[j],]
                               
                               # sigma_spat doesn't matter prior to eigen-decomposition
                               image_cor_mat <- get_correlation(
                                   fov_data = image_dat,
                                   zero_distance = zero_dist_use,
                                   sigma_spat = 1,
                                   cor_type = 'sqexp',
                                   induce_sparsity = FALSE
                               )
                               
                               eigens <- eigen(image_cor_mat)
                               
                               # weight eigenvectors by proportion of variance explained
                               offset <- colSums(t(eigens$vectors) * (eigens$values / sum(eigens$values)))
                               # standardize
                               (offset - mean(offset)) / sd(offset)
                               
                           }

parallel::stopCluster(cl = my_cluster)


sim_data$offset <- spatial_offsets
sim_data<- sim_data[order(sim_data$outcome / sim_data$numNeighbors),]

p1 <- ggplot(sim_data, aes(x = x, y = y, col = outcome/numNeighbors)) + 
    geom_point(size = 1) + 
    scale_color_gradient(low = 'grey80', high = 'red') + 
    facet_grid(~subjectID) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.text = element_text(size = 12)) +
    labs(color = "Outcome") + 
    scale_x_continuous(breaks = seq(0, 750, 250))+ 
    scale_y_continuous(breaks = seq(0, 750, 250))

p2 <- ggplot(sim_data, aes(x = x, y = y, col = offset)) + 
    geom_point(size = 1) + 
    scale_color_viridis_c() + 
    facet_grid(~subjectID) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.text = element_text(size = 12)) +
    labs(color = "Spatial    \noffset") + 
    scale_x_continuous(breaks = seq(0, 750, 250))+ 
    scale_y_continuous(breaks = seq(0, 750, 250))

png('figures/sim_datasets_offset.png', units = 'in', res = 500, height =6, width = 10)
grid.arrange(p1, p2)
dev.off()

################################################################################
### Figure 1 - example datasets

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
# Figure 2 - bias and coverage of 95% CI

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
                                         'no_corr',
                                         'pc_sqexp'),
                              labels = c('Full model',
                                         'No spatial correlation',
                                         'Eigen-decomposition'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))


n_subjects_plot <- 30
n_image_sub_plot <- 5
zero_distance_plot <- 100


# Bias 
p1 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot & 
                        zero_distance == zero_distance_plot), 
             aes(x = exp(beta_val), y = bias_rel,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'tomato', 'dodgerblue')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Relative Bias',
         col = 'Model') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1) + 
    ylim(-0.4, 0.4)


# coverage
p2 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot & 
                        zero_distance == zero_distance_plot), 
             aes(x = exp(beta_val), y = cover,  col = model_type)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'tomato', 'dodgerblue')) + 
    ggtitle('Coverage of 95% CIs') +
    labs(x = 'Odds ratio for responders vs non-responders',
         y = 'Coverage  probability',
         col = 'Model') +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1)

# Extract the legend from one plot
legend <- get_legend(
    p1 + theme(legend.position = "right")
)


# Remove legends from the individual panels
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")

# Combine plots vertically, add shared legend to the right


png('figures/bias_cover.png', units = 'in', res = 500, height =7, width = 12)
plot_grid(
    plot_grid(p1_noleg, p2_noleg, ncol = 1, align = "v"),
    legend,
    rel_widths = c(1, 0.3)
)
dev.off()




################################################################################
# Figure 3 - power by dataset size

or_stats$n_image_sub <- factor(or_stats$n_image_sub, 
                               labels = c('1 image per subject',
                                          '5 images per subject'))
or_stats <- or_stats[order(or_stats$beta_idx, or_stats$n_subjects, or_stats$zero_distance,
                           or_stats$sigma_spat, or_stats$model_type, or_stats$n_image_sub),]

or_stats$coef <- factor(or_stats$coef,
                        levels = c("OR_R", 'OR_FOV'),
                        labels = c('Subject-level predictor',
                                   'Image-level predictor'))
or_stats$n_subjects <- as.numeric(as.character(or_stats$n_subjects))

png('figures/power_datasize.png', units = 'in', res = 500, height =3.5, width = 9)
ggplot(subset(or_stats, 
              beta_idx == 7  & 
                  zero_distance == 100 &
                  sigma_spat_fac == 'High spatial correlation' &
                  model_type %in% c('Eigen-decomposition')), 
       aes(x = n_subjects, y = power, col = n_image_sub)) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    geom_line(linewidth = 1.5) +
    facet_nested(  ~ coef) +
    labs(x = 'Number of subjects',
         y = 'Power',
         col = '') +
    ylim(0, 1) +
    theme(legend.key.width = unit(2,"cm")) + 
    scale_color_manual(values = c('#1B9E77',
                                  "#D95F02"),
                       guide = guide_legend(reverse = TRUE) ) 
dev.off()


# png('figures/power_datasize.png', units = 'in', res = 500, height =6, width = 12)
# ggplot(subset(or_stats, 
#               beta_idx > 5 & 
#                   zero_distance == 100 &
#                   model_type %in% c('Eigen-decomposition')), 
#        aes(x = exp(beta_val), y = power, col = n_image_sub,  linetype = n_subjects)) +
#     # geom_point(size = 2) +
#     geom_line(linewidth = 1) +
#     facet_nested(coef  ~ sigma_spat_fac) +
#     labs(x = 'Odds ratio for subject-level predictor',
#          y = 'Power',
#          linetype = 'Number of subjects',
#          col = '') +
#     ylim(0, 1) +
#     theme(legend.key.width = unit(2,"cm")) + 
#     scale_color_manual(values = c('#1B9E77',
#                                   "#D95F02"))
# dev.off()

################################################################################
# Table 1 - MSE by dataset size
# time to compute

time_tab <- readRDS('results/time_tab.rds')
ncells_summary <- readRDS('results/ncells_summary.rds')

time_tab$median_IQR <- paste0(sprintf("%.2f", round(time_tab$median_time/60, 2)), 
                              ' (',
                              sprintf("%.2f", round(time_tab$lower_time/60, 2)),
                              ', ',
                              sprintf("%.2f", round(time_tab$upper_time/60, 2)),
                              ')')


time_tab <- time_tab[,c('n_subjects', 'n_image_sub',
                        'model_type', 'median_IQR')]

time_tab_wide <- pivot_wider(time_tab,
                             id_cols = c(n_subjects, n_image_sub),
                             names_from=model_type,
                             values_from=median_IQR) %>%
    data.frame()
time_tab_wide$median_cells <- ncells_summary$median_cells
time_tab_wide <- time_tab_wide[,c('n_subjects', 'n_image_sub', 'median_cells',
                                  'inla','no_corr',  'pc_sqexp')]

colnames(time_tab_wide)[4:6] <- c('Full model',
                                  'No spatial correlation',
                                  'Eigen-decomposition')

time_tab_wide$median_cells <- scales::comma(time_tab_wide$median_cells)

time_tab_wide

kable(time_tab_wide) %>%
    collapse_rows(columns = 1, valign = "middle")


