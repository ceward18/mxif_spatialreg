
fit_model <- function(model_type, sim_data, model_info, beta_val, beta2) {
    
    zero_distance <- model_info$zero_distance
    
    if (model_type == 'no_corr') {
        # MODEL with spatial correlation ignored
        
        start_model <- Sys.time()
        
        # remove those with no neighbors
        sub_dat <- sim_data[sim_data$numNeighbors > 0,]
        
        spfit <- glmmTMB(cbind(outcome, numNeighbors-outcome) ~ 
                             response + fov_type +
                             (1|subjectID) +
                             (1|imageID),
                         data = sub_dat,
                         family = binomial,
                         REML = TRUE)
        
        end_model <- Sys.time()
        
        model_fit_time <-  as.numeric(end_model - start_model, units = 'secs')
        
        get_summary_table(spfit, sim_data, beta_val, beta2,
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
        
        sigma_spat_use <- sqrt(mean(vario_ests[,1]))
        zero_dist_use <- mean(vario_ests[,2])
        
        # print(sigma_spat_use)
        
        # get correlation type to use - exponential or squared exponential
        cor_type_use <- gsub('.*_', '', model_type)
        # print(cor_type_use)
        
        
        spatial_offset <- foreach(j = 1:length(unique_images), 
                                  .combine = 'c') %dopar% {
                                      
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
                                      # standardize
                                      (offset - mean(offset)) / sd(offset)
                                      
                                  }
        
        
        parallel::stopCluster(cl = my_cluster)
        
        end_eigen <- Sys.time()
        
        eigen_decomp_time <-  as.numeric(end_eigen - start_eigen, units = 'secs')
        
        sim_data$spatial_offset <- sigma_spat_use * spatial_offset
        
        # PC model fitting
        start_model <- Sys.time()
        
        # remove those with no neighbors
        sub_dat <- sim_data[sim_data$numNeighbors > 0,]
        
        spfit <- glmmTMB(cbind(outcome, numNeighbors-outcome) ~ 
                             response + fov_type + 
                             (1|subjectID) +
                             (1|imageID) + 
                             offset(spatial_offset),
                         data = sub_dat,
                         family = binomial,
                         REML = TRUE)
        
        
        end_model <- Sys.time()
        
        model_fit_time <- as.numeric(end_model - start_model, units = 'secs')
        
        
        get_summary_table(spfit, sub_dat, beta_val, beta2,
                          sigma_spat_use, model_type = model_type,
                          model_info,
                          model_fit_time, eigen_decomp_time)
        
    }
    
    
    
}
