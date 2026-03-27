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
# Figure 1 - example datasets and contribution of spatial correlation from normalized offset

rho_25 <- optim(0, function(x) (exp(-25^2 / (2 * x^2)) - 0.001)^2, 
                method = 'Brent',
                lower = 0, upper = 500)$par

rho_100 <- optim(0, function(x) (exp(-100^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par

rho_250 <- optim(0, function(x) (exp(-250^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par

set.seed(100)
sim_data1 <- sim_data_fn(betas = c( logit(0.10), 0, 0), 
                         nSub = 2, 
                         nImagePerSub = 2, 
                         rho = rho_25,
                         sigma_spat = 2,    # spatial SD
                         sigma_sub = 0.4,            # between subjects SD
                         sigma_image = 0.4)          # between images SD
sim_data2 <- sim_data_fn(betas = c( logit(0.10), 0, 0), 
                         nSub = 2, 
                         nImagePerSub = 2, 
                         rho = rho_100,
                         sigma_spat = 2,    # spatial SD
                         sigma_sub = 0.4,            # between subjects SD
                         sigma_image = 0.4)          # between images SD
sim_data3 <- sim_data_fn(betas = c( logit(0.10), 0, 0), 
                         nSub = 2, 
                         nImagePerSub = 2, 
                         rho = rho_250,
                         sigma_spat = 2,    # spatial SD
                         sigma_sub = 0.4,            # between subjects SD
                         sigma_image = 0.4)          # between images SD

get_eigen_offset <- function(sim_data) {
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
    sim_data
}

sim_data1 <- get_eigen_offset(sim_data1)
sim_data2 <- get_eigen_offset(sim_data2)
sim_data3 <- get_eigen_offset(sim_data3)

# just keep one image from each
sim_data <- rbind.data.frame(subset(sim_data1, imageID == 1),
                             subset(sim_data2, imageID == 2),
                             subset(sim_data3, imageID == 3))

sim_data<- sim_data[order(sim_data$outcome / sim_data$numNeighbors),]

sim_data$imageID <- factor(sim_data$imageID,
                           labels = c('Short range', 
                                      'Medium range', 
                                      'Long range'))

fig1_theme <- theme_bw() + 
    theme(strip.text = element_text(size = 20),
          strip.background = element_rect(fill = 'white'),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank())

p1 <- ggplot(sim_data, aes(x = x, y = y, col = outcome/numNeighbors)) + 
    geom_point(size = 2) + 
    scale_color_gradient(low = 'grey80', high = 'red') + 
    facet_grid(~imageID) +
    fig1_theme +
    guides(col = guide_colorbar(barheight = unit(6, "cm"))) +
    labs(color = "Outcome\n", x = '', y = '') 

p2 <- ggplot(sim_data, aes(x = x, y = y, col = offset)) + 
    geom_point(size = 2) + 
    scale_color_viridis_c() + 
    facet_grid(~imageID) +
    fig1_theme +
    guides(col = guide_colorbar(barheight = unit(6, "cm"))) +
    labs(color = expression(tilde(S)[ijk]~'         '), x = '', y = '')  

fig1 <- grid.arrange(p1, p2)


# Save the plot as an EPS file
ggsave("figures/fig1.eps", plot = fig1, device = "eps", 
       width = 12, height = 8, units = "in")


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
                                         'pc_sqexp', 
                                         'no_corr'),
                              labels = c('INLA-SPDE',
                                         'Eigen-decomposition',
                                         'No spatial correlation'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))


n_subjects_plot <- 30
n_image_sub_plot <- 5
zero_distance_plot <- 100



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


# Bias 
p1 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot & 
                        zero_distance == zero_distance_plot), 
             aes(x = exp(beta_val), y = bias_rel,  col = model_type)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1) + 
    # geom_point(size = 3, alpha = 0.8) + 
    geom_line(linewidth = 1.5, alpha = 0.8) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('goldenrod2', 'dodgerblue', 'tomato')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Relative Bias',
         col = 'Model') +
    fig2_theme + 
    scale_y_continuous(limits = c(-0.4, 0.4))

# coverage
p2 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot & 
                        zero_distance == zero_distance_plot), 
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

# Combine plots vertically, add shared legend to the right
fig2 <- plot_grid(
    plot_grid(p1_noleg, p2_noleg, ncol = 1, align = "v"),
    legend,
    rel_widths = c(1, 0.3)
)

# Save the plot as an EPS file
ggsave("figures/fig2.eps", plot = fig2, device = "pdf", 
       width = 12, height = 8, units = "in")




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

fig3 <- ggplot(subset(or_stats, 
                      beta_idx == 7  & 
                          zero_distance == 100 &
                          sigma_spat_fac == 'High spatial correlation' &
                          model_type %in% c('Eigen-decomposition')), 
               aes(x = n_subjects, y = power, col = n_image_sub)) +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 1) +
    geom_line(linewidth = 1.5) +
    facet_nested(  ~ coef) +
    labs(x = 'Number of subjects',
         y = 'Power',
         col = '') +
    ylim(0, 1) +
    theme_bw() + 
    theme(strip.text = element_text(size = 16),
          strip.background = element_rect(fill = 'white'),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12),
          legend.key.width = unit(1.5,"cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + 
    scale_color_manual(values = c('#1B9E77',
                                  "#D95F02"),
                       guide = guide_legend(reverse = TRUE) ) 

ggsave("figures/fig3.eps", plot = fig3, device = "eps", 
       width = 9, height = 3.5, units = "in")



################################################################################
# Figure 4 - Compute time by dataset size

time_n_cells <- readRDS('results/time_n_cells.rds')


time_n_cells$model_type <- factor(time_n_cells$model_type,
                                  levels = c('inla',
                                             'pc_sqexp', 
                                             'no_corr'),
                                  labels = c('INLA-SPDE',
                                             'Eigen-decomposition',
                                             'No spatial correlation'))

fig4 <- ggplot(time_n_cells, aes(x = median_cells, y = median_time/60, col = model_type)) + 
    geom_line() + 
    geom_point() + 
    labs(y = 'Median compute time (min)', 
         x = 'Number of index cells',
         col = 'Model')+
    scale_x_continuous(labels = scales::comma, breaks = seq(0, 300000, 50000))+
    scale_color_manual(values = c('goldenrod2', 'dodgerblue', 'tomato')) +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())


# Save the plot as an EPS file
ggsave("figures/fig4.eps", plot = fig4, device = "eps", 
       width = 6, height = 3, units = "in")


# linear rate of change (for discussion and abstract)
inla_mod <- lm(median_time ~ mean_cells, data = subset(time_n_cells, model_type == 'INLA-SPDE'))
inla_mod$coefficients[2] * 50000 / 60 # in minutes


eigen_mod <- lm(median_time ~ mean_cells, data = subset(time_n_cells, model_type == 'Eigen-decomposition'))
eigen_mod$coefficients[2] * 50000  # in seconds


################################################################################
# Supplemental figure 1 - variability of S as a function of n

rho_100 <- optim(0, function(x) (exp(-100^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par

set.seed(100)
sim_data1 <- sim_data_fn(betas = c( logit(0.10), 0, 0), 
                         nSub = 10, 
                         nImagePerSub = 1, 
                         rho = rho_100,
                         sigma_spat = 2,    # spatial SD
                         sigma_sub = 0.4,            # between subjects SD
                         sigma_image = 0.4)          # between images SD


get_eigen_offset_unstandardized <- function(sim_data) {
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
                                   # standardize (want unstandardized just for illustration!)
                                   # (offset - mean(offset)) / sd(offset)
                                   
                               }
    
    parallel::stopCluster(cl = my_cluster)
    
    
    sim_data$offset <- spatial_offsets
    sim_data
}

sim_data1 <- get_eigen_offset_unstandardized(sim_data1)


# number of cells per image
# images to use for illustration of the sample size effect: 7, 6, 9
n_cells_image <- sim_data1 %>%
    group_by(imageID) %>%
    summarize(n_cells_image = length(tumorID))
sim_data1 <- merge(sim_data1, n_cells_image, by  = 'imageID', all.x = T)

sim_data1$n_cells_image <- paste0('K = ', sim_data1$n_cells_image)
sim_data1$n_cells_image <- factor(sim_data1$n_cells_image,
                                  levels = paste0('K = ', sort(n_cells_image$n_cells_image)))


png('figures/fig_s1_s_std.png', units = 'in', res = 500, height =3, width = 7)
ggplot(subset(sim_data1, imageID %in% c(7, 6, 9)),
       aes(x = offset)) + 
    geom_histogram(bins = 40) + 
    facet_wrap(~n_cells_image) + 
    theme_bw() +
    theme(strip.text = element_text(size = 12),
          strip.background = element_rect(fill = 'white'),
          axis.title = element_text(size = 12), 
          axis.text = element_text(size = 8)) + 
    labs(x = expression(S[ijk]), y = 'Frequency') + 
    coord_cartesian(xlim = c(-0.015, 0.015))
dev.off()




################################################################################
# Supplemental figure 2 - exponential versus squared exponential

or_stats <- subset(prop_stats, 
                   model_type %in% c('pc_exp', 
                                     'pc_sqexp'))


or_stats$model_type <- factor(or_stats$model_type,
                              levels = c('pc_exp',
                                         'pc_sqexp'),
                              labels = c('Exponential\n(misspecified)',
                                         'Squared exponential\n(correct)'))

or_stats$sigma_spat_fac <- factor(or_stats$sigma_spat,
                                  labels = c('Low spatial correlation',
                                             'Medium spatial correlation',
                                             'High spatial correlation'))


or_stats$zero_distance <- factor(or_stats$zero_distance,
                                 labels = c('Shorter range',
                                            'Longer range'))



n_subjects_plot <- 30
n_image_sub_plot <- 5
zero_distance_plot <- 100


figs1_theme <- theme_bw() + 
    theme(strip.text = element_text(size = 16),
          strip.background = element_rect(fill = 'white'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.width = unit(1, "cm"),   
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.title = element_text(size = 16), 
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 18, h = 0.5))


# Bias 
p1 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot ), 
             aes(x = exp(beta_val), y = bias_rel,  col = model_type, linetype = zero_distance)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue')) + 
    ggtitle('Bias') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Relative Bias',
         col = 'Model',
         linetype = 'Range of spatial correlation') +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 1)  + 
    figs1_theme


# coverage
p2 <- ggplot(subset(or_stats, coef == 'OR_R' & 
                        n_subjects == n_subjects_plot & 
                        n_image_sub == n_image_sub_plot), 
             aes(x = exp(beta_val), y = cover,  col = model_type, linetype = zero_distance)) +
    geom_point(size = 2) + 
    geom_line(linewidth = 1) +
    facet_nested( ~ sigma_spat_fac ) +
    scale_color_manual(values = c('tomato', 'dodgerblue')) + 
    ggtitle('Coverage of 95% CIs') +
    labs(x = 'Odds ratio for subject-level predictor',
         y = 'Coverage  probability',
         col = 'Model',
         linetype = 'Range of spatial correlation') +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1) + 
    figs1_theme

# Extract the legend from one plot
legend <- get_legend(
    p1 + theme(legend.position = "right")
)


# Remove legends from the individual panels
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")

# Combine plots vertically, add shared legend to the right


png('figures/fig_s2_bias_cover_cov.png', units = 'in', res = 500, height =7, width = 12)
plot_grid(
    plot_grid(p1_noleg, p2_noleg, ncol = 1, align = "v"),
    legend,
    rel_widths = c(1, 0.3)
)
dev.off()
