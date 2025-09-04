
fit_model <- function(model_type, sim_data, model_info, beta0, beta_val) {
    
    zero_distance <- model_info$zero_distance
    
    if (model_type == 'no_corr') {
        # MODEL with spatial correlation ignored
        
        start <- Sys.time()
        
        # remove those with no neighbors
        sub_dat <- sim_data[sim_data$numNeighbors > 0,]
        
        spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ response +
                           (1|subjectID) +
                           (1|imageID),
                       data = sub_dat,
                       family = binomial,
                       method = 'REML')
        
        end <- Sys.time()
        
        model_fit_time <-  as.numeric(end - start, units = 'secs')
        
        get_summary_table(spfit, sim_data, beta0, beta_val,
                          sigma_spat_use, model_type,
                          model_info,
                          model_fit_time, 
                          eigen_decomp_time = 0)
        
        
    } else {
        # principal components models
        
        start_eigen <- Sys.time()
        
        unique_images <- unique(sim_data$imageID)
        
        # estimate variograms for each image in parallel
        my_cluster <- makeCluster(10)
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
        
        # print(sigma_spat_use)
        
        # get correlation type to use - exponential or squared exponential
        cor_type_use <- gsub('.*_', '', model_type)
        # print(cor_type_use)
        
        
        spatial_offsets <- foreach(j = 1:length(unique_images), 
                                   .combine = 'rbind') %dopar% {
                                       
                                       source('helper_functions.R')
                                       
                                       image_dat <- sim_data[sim_data$imageID == unique_images[j],]
                                       
                                       # sigma_spat doesn't matter prior to eigen-decomposition
                                       image_cor_mat <- get_correlation(
                                           fov_data = image_dat,
                                           zero_distance = zero_dist_use,
                                           sigma_spat = 1,
                                           cor_type = cor_type_use,
                                           induce_sparsity = FALSE
                                       )
                                       
                                       eigens <- eigen(image_cor_mat)
                                       
                                       # weight eigenvectors by proportion of variance explained
                                       offset <- colSums(t(eigens$vectors) * (eigens$values / sum(eigens$values)))
                                       # standardize?
                                       norm_offset <- (offset - mean(offset)) / sd(offset)
                                       
                                       # multiply by sigma estimate individually
                                       sigma_offset <- offset * vario_ests[j,1]
                                       
                                       
                                       # multiply by sigma estimate individually
                                       sigma_norm_offset <- norm_offset * vario_ests[j,1]
                                       
                                       
                                       cbind(offset,
                                             sigma_offset,
                                             norm_offset,
                                             sigma_norm_offset)
                                       
                                   }
        
        
        parallel::stopCluster(cl = my_cluster)
        
        end_eigen <- Sys.time()
        
        sim_data$offset <- spatial_offsets[,1]
        sim_data$sigma_offset <- spatial_offsets[,2]
        sim_data$norm_offset <- spatial_offsets[,3]
        sim_data$sigma_norm_offset <- spatial_offsets[,4]
        
        eigen_decomp_time <-  as.numeric(end_eigen - start_eigen, units = 'secs')
        
        
        # PC models fitting
        pc_model_combos <- expand.grid(offset_base = c('raw', 'norm'),
                    sigma_adjust = c('image', 'overall'),
                    sigma_type = c('1', 'sigma'))
        
        # if no adjustment doesn't matter if done within image or not
        pc_model_combos <- pc_model_combos[-which(pc_model_combos$sigma_type == '1' & 
                                               pc_model_combos$sigma_adjust == 'image'),]
        pc_model_combos <- pc_model_combos[order(pc_model_combos$offset_base, 
                                                 pc_model_combos$sigma_adjust, 
                                                 pc_model_combos$sigma_type),]
        
        pc_model_types <- paste0(model_type, '_', pc_model_combos$offset_base, '_', 
                                 pc_model_combos$sigma_adjust, '_', pc_model_combos$sigma_type)
        
        sum_tabs_list <- list()
        
        for (k in 1:length(pc_model_types)) {
            
            print(pc_model_types[k])
            
            start <- Sys.time()
            
            # remove those with no neighbors
            sub_dat <- sim_data[sim_data$numNeighbors > 0,]
            
            # configure spatial offset according to model type
            
            # use raw or normalized eigenvectors
            if (grepl('raw', pc_model_types[k])) {
                
                if (grepl('image', pc_model_types[k])) {
                    # within image sigma adjustments
                    sub_dat$spatial_offset <- sub_dat$sigma_offset 
                    
                } else {
                    # overall adjustments
                    
                    if (grepl('sigma$', pc_model_types[k])) {
                        sub_dat$spatial_offset <- sigma_spat_use * sub_dat$offset 
                    } else {
                        sub_dat$spatial_offset <- sub_dat$offset 
                        
                    }
                }
            } else if (grepl('norm', pc_model_types[k])) {
                
                if (grepl('image', pc_model_types[k])) {
                    # within image sigma adjustments
                    sub_dat$spatial_offset <- sub_dat$sigma_norm_offset 
                    
                } else {
                    # overall adjustments
                    
                    if (grepl('sigma$', pc_model_types[k])) {
                        sub_dat$spatial_offset <- sigma_spat_use * sub_dat$norm_offset 
                    } else {
                        sub_dat$spatial_offset <- sub_dat$norm_offset 
                        
                    }
                }
            }
            
            
            spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ response +
                               (1|subjectID) +
                               (1|imageID) + 
                               offset(spatial_offset),
                           data = sub_dat,
                           family = binomial,
                           method = 'REML')
            
            
            
            end <- Sys.time()
            
            model_fit_time <- as.numeric(end - start, units = 'secs')
            
            
            sum_tabs_list[[k]] <- get_summary_table(spfit, sub_dat, beta0, beta_val,
                                                    sigma_spat_use, model_type = pc_model_types[k],
                                                    model_info,
                                                    model_fit_time, eigen_decomp_time)
            
        }
        
        do.call('rbind.data.frame', sum_tabs_list)
        
    }
    
    
    
}
