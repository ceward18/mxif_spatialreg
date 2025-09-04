################################################################################
# combine batch runs
################################################################################

outputFolder <- 'output'
resultsFolder <- 'results'

outputFiles <- sort(list.files(paste0('./', outputFolder)))

# Three batches 
# for each batch convert to summary stats then combine all summary stats

res_all <- readRDS(paste0('./', outputFolder, '/', outputFiles[1]))

for (i in 2:length(outputFiles)) {
    if (i %% 25 == 0) print(i)
    res_i <- readRDS(paste0('./', outputFolder, '/', outputFiles[i]))
    res_all <-rbind.data.frame(res_all, res_i)
}

saveRDS(res_all, paste0('./', resultsFolder, '/res_all.rds'))


library(dplyr)

# MSE, absolute difference, coverage probability, avg time
prop_stats <- res_all %>%
    group_by(zero_distance, sigma_spat, coef, model_type, beta_idx) %>%
    summarise(nsims = length(sim_number),
              mse = mean((est - truth)^2),
              mse_rel = mean((est/truth - truth/truth)^2),
              bias = mean((est - truth)),
              bias_rel = mean((est - truth)/truth),
              abs_diff = mean(abs(est - truth)),
              abs_diff_rel = mean(abs(est/truth - truth/truth)),
              cover = mean(truth > lower & truth < upper),
              time = mean(time),
              truth = mean(truth)) %>%
    data.frame()


saveRDS(prop_stats, paste0('./', resultsFolder, '/prop_stats.rds'))


time_tab <- subset(res_all, coef == 'OR' & 
                       model_type %in% c('no_corr', 
                                         'pc_sqexp_norm_sigma',
                                         'nngp')) %>%
    group_by(model_type) %>%
    summarise(avg_time = mean(time),
              median_time = median(time),
              lower_time = quantile(time, 0.25),
              upper_time = quantile(time, 0.75)) %>%
    data.frame()

saveRDS(time_tab, paste0('./', resultsFolder, '/time_tab.rds'))



# median number of tumor cells by n_subjects and n_images
ncells_summary <- subset(res_all, coef == 'OR' & 
           model_type %in% c('no_corr')) %>%
    summarise(mean_cells = mean(n_cells),
              median_cells = round(median(n_cells)),
              lower_cells = quantile(n_cells, 0.25),
              upper_cells = quantile(n_cells, 0.75)) %>%
    data.frame()

saveRDS(ncells_summary, paste0('./', resultsFolder, '/ncells_summary.rds'))
