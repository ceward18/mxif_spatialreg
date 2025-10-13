################################################################################
# combine batch runs
################################################################################

outputFolder <- 'output'
resultsFolder <- 'results'

outputFiles <- sort(list.files(paste0('./', outputFolder)))
outputFiles <- outputFiles[grep('res_batch', outputFiles)]

# any missing?
all_files <- paste0('res_batch_', sprintf("%04d",1:1620), '.rds')
all_files[!all_files %in% outputFiles]

# Three batches 
# for each batch convert to summary stats then combine all summary stats


outputFiles1 <- outputFiles[1:540]

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles1[1]))

for (i in 2:length(outputFiles1)) {
    if (i %% 25 == 0) print(i)
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles1[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}

saveRDS(res_all, paste0('./', resultsFolder, '/res_all_1.rds'))


# second batch
outputFiles2 <- outputFiles[541:1080]

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles2[1]))

for (i in 2:length(outputFiles2)) {
    if (i %% 25 == 0) print(i)
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles2[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}

saveRDS(res_all, paste0('./', resultsFolder, '/res_all_2.rds'))




# third batch
outputFiles3 <- outputFiles[1081:length(outputFiles)]

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles3[1]))

for (i in 2:length(outputFiles3)) {
    if (i %% 25 == 0) print(i)
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles3[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}

saveRDS(res_all, paste0('./', resultsFolder, '/res_all_3.rds'))

# inla batch
outputFiles <- sort(list.files(paste0('./', outputFolder)))
outputFiles <- outputFiles[grep('res_inla', outputFiles)]

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles[1]))

for (i in 2:length(outputFiles)) {
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}
saveRDS(res_all, paste0('./', resultsFolder, '/res_all_inla.rds'))


################################################################################
### combine all and summarize

resultsFolder <- 'results'

res_all_1 <- readRDS(paste0('./', resultsFolder, '/res_all_1.rds'))
res_all_2 <- readRDS(paste0('./', resultsFolder, '/res_all_2.rds'))
res_all_3 <- readRDS(paste0('./', resultsFolder, '/res_all_3.rds'))
res_all_inla <- readRDS(paste0('./', resultsFolder, '/res_all_inla.rds'))
res_all_inla$zero_distance <- 100


res_all <- rbind.data.frame(res_all_1, res_all_2, res_all_3, res_all_inla)

# forgot to exponentiate true value
res_all$truth[which(res_all$coef == 'OR_FOV' & res_all$model_type != 'inla')] <-
    exp(res_all$truth[which(res_all$coef == 'OR_FOV' & res_all$model_type != 'inla')])


library(dplyr)

# MSE, absolute difference, coverage probability, avg time
prop_stats <- res_all %>%
    group_by(n_subjects, n_image_sub, zero_distance, sigma_spat, 
             coef, model_type, beta_idx) %>%
    summarise(nsims = length(sim_number),
              mse = mean((est - truth)^2),
              mse_rel = mean((est/truth - truth/truth)^2),
              bias = mean((est - truth)),
              bias_rel = mean((est - truth)/truth),
              abs_diff = mean(abs(est - truth)),
              abs_diff_rel = mean(abs(est/truth - truth/truth)),
              cover = mean(truth > lower & truth < upper, na.rm = T),
              power = mean((lower > 1) | (upper < 1), na.rm = T), 
              time = mean(time),
              truth = mean(truth)) %>%
    data.frame()


saveRDS(prop_stats, paste0('./', resultsFolder, '/prop_stats.rds'))



time_tab <- subset(res_all, coef == 'OR_R' & 
                       model_type %in% c('no_corr', 
                                         'pc_sqexp', 
                                         'inla')) %>%
    group_by(model_type, n_subjects, n_image_sub) %>%
    summarise(avg_time = mean(time),
              median_time = median(time),
              lower_time = quantile(time, 0.25),
              upper_time = quantile(time, 0.75),
              median_eigen_time = median(eigen_decomp_time),
              median_model_time = median(model_fit_time),
              median_summary_time = median(summary_time)) %>%
    data.frame()

saveRDS(time_tab, paste0('./', resultsFolder, '/time_tab.rds'))



# median number of tumor cells by n_subjects and n_images
ncells_summary <- subset(res_all, coef == 'OR_R' & 
           model_type %in% c('no_corr')) %>%
    group_by(n_subjects, n_image_sub) %>%
    summarise(mean_cells = mean(n_cells),
              median_cells = round(median(n_cells)),
              lower_cells = quantile(n_cells, 0.25),
              upper_cells = quantile(n_cells, 0.75)) %>%
    data.frame()

saveRDS(ncells_summary, paste0('./', resultsFolder, '/ncells_summary.rds'))
