
fit_model <- function(model_type, sim_data, model_info, beta0, beta_val) {
    
    zero_distance <- model_info$zero_distance
    
    if (model_type == 'no_corr') {
        # MODEL with spatial correlation ignored
        
        start <- Sys.time()
        
        # remove those with no neighbors
        sub_dat <- sim_data[sim_data$numNeighbors > 0,]
        
        spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ marker,
                       data = sub_dat,
                       family = binomial,
                       method = 'REML')
        
        end <- Sys.time()
        
        model_fit_time <-  as.numeric(end - start, units = 'secs')
        
        get_summary_table(spfit, sim_data, beta0, beta_val,
                          sigma_spat_use=NA, zero_dist_use=NA,
                          model_type,
                          model_info,
                          model_fit_time, 
                          eigen_decomp_time = 0)
        
        
    } else if (grepl('pc', model_type)) {
        # principal components models
        
        start_eigen <- Sys.time()
        
        # estimate variogram
        vario_ests <-  variogram_est(fov_data = sim_data)
        sigma_spat_use <- mean(vario_ests[1])
        zero_dist_use <- mean(vario_ests[2])
        
        # get correlation type to use - exponential or squared exponential
        cor_type_use <- gsub('.*_', '', model_type)
        # print(cor_type_use)
        
        # sigma_spat doesn't matter prior to eigen-decomposition
        image_cor_mat <- get_correlation(
            fov_data = sim_data,
            zero_distance = zero_dist_use,
            cor_type = cor_type_use
        )
        
        eigens <- eigen(image_cor_mat)
        
        # weight eigenvectors by proportion of variance explained
        offset <- colSums(t(eigens$vectors) * (eigens$values / sum(eigens$values)))
        # normalize?
        norm_offset <- (offset - mean(offset)) / sd(offset)
        
        # multiply by sigma estimate 
        sigma_offset <- offset * sigma_spat_use
        
        # multiply by sigma estimate
        sigma_norm_offset <- norm_offset * sigma_spat_use
        
        sim_data$offset <- offset
        sim_data$sigma_offset <- sigma_offset
        sim_data$norm_offset <- norm_offset
        sim_data$sigma_norm_offset <- sigma_norm_offset
        
        end_eigen <- Sys.time()
        
        eigen_decomp_time <-  as.numeric(end_eigen - start_eigen, units = 'secs')
        
        
        # PC models fitting
        pc_model_combos <- expand.grid(offset_base = c('raw', 'norm'),
                                       sigma_type = c('1', 'sigma'))
        
        # if no adjustment doesn't matter if done within image or not
        pc_model_combos <- pc_model_combos[order(pc_model_combos$offset_base, 
                                                 pc_model_combos$sigma_type),]
        
        pc_model_types <- paste0(model_type, '_', pc_model_combos$offset_base, '_', 
                                 pc_model_combos$sigma_type)
        
        sum_tabs_list <- list()
        
        for (k in 1:length(pc_model_types)) {
            
            print(pc_model_types[k])
            
            start <- Sys.time()
            
            # remove those with no neighbors
            sub_dat <- sim_data[sim_data$numNeighbors > 0,]
            
            # configure spatial offset according to model type
            
            # use raw or normalized eigenvectors
            if (grepl('raw', pc_model_types[k])) {
                
                if (grepl('sigma$', pc_model_types[k])) {
                    sub_dat$spatial_offset <- sub_dat$sigma_offset 
                } else {
                    sub_dat$spatial_offset <- sub_dat$offset 
                }
                
            } else if (grepl('norm', pc_model_types[k])) {
                
                if (grepl('sigma$', pc_model_types[k])) {
                    sub_dat$spatial_offset <- sub_dat$sigma_norm_offset 
                } else {
                    sub_dat$spatial_offset <- sub_dat$norm_offset   
                }
            }
            
            
            spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ marker +
                               offset(spatial_offset),
                           data = sub_dat,
                           family = binomial,
                           method = 'REML')
            
            end <- Sys.time()
            
            model_fit_time <- as.numeric(end - start, units = 'secs')
            
            
            sum_tabs_list[[k]] <- get_summary_table(spfit, sub_dat, beta0, beta_val,
                                                    sigma_spat_use, zero_dist_use,
                                                    model_type = pc_model_types[k],
                                                    model_info,
                                                    model_fit_time, eigen_decomp_time)
            
        }
        
        do.call('rbind.data.frame', sum_tabs_list)
        
    } else {
        # NNGP model
        
        start_nngp <- Sys.time()
        
        # convert to probability and use normal distribution because NNGP can't
        # do binomial when size > 1
        sim_data$prop <- sim_data$outcome / sim_data$numNeighbors
        sim_data$logit_prop <- log(sim_data$prop / (1 - sim_data$prop)) 
        sim_data$logit_prop[sim_data$prop == 0] <- log(0.01 / 0.99) 
        sim_data$logit_prop[sim_data$prop == 1] <- log(0.99 / 0.01) 
        # remove those with no neighbors
        sub_dat <- sim_data[sim_data$numNeighbors > 0,]
        
        # use variogram for initial values
        vario_ests <-  variogram_est(fov_data = sim_data)
        sigma_spat_use <- mean(vario_ests[1])
        zero_dist_use <- mean(vario_ests[2])
        
        
        # initial values
        starting <- list("phi" = sqrt(-log(0.001))/zero_dist_use,
                         "sigma.sq" = sigma_spat_use,
                         "tau.sq" = 0.01)
        
        # proposal standard deviations
        tuning <- list("phi" = 0.1,
                       "sigma.sq" = 0.05,
                       "tau.sq" = 0.02)
        
        # priors
        priors <- list("phi.Unif" = c(sqrt(-log(0.001))/200, sqrt(-log(0.001))/5),
                       "sigma.sq.IG" = c(2, sigma_spat_use),
                       "tau.sq.IG" = c(2, 0.01))
        
        # covariance function
        cov.model <- "gaussian"
        
        
        # fit the model
        nngp_mod <- spNNGP(logit_prop ~ marker ,
                           coords = c("x", "y"),
                           data = sub_dat,
                           n.neighbors = 10, 
                           starting = starting,
                           tuning = tuning,
                           priors = priors,
                           cov.model = cov.model,
                           n.samples = 20000,                # number of iterations
                           n.report = 5000,                  # progress updates
                           fit.rep = TRUE,
                           sub.sample = list(start = 5000), # burn-in
                           n.omp.threads = 10)
        
        
        end_nngp <- Sys.time()
        
        model_fit_time <- as.numeric(end_nngp - start_nngp, units = 'secs')
        
        
        get_summary_table_nngp(nngp_mod, beta0, beta_val,
                               model_type,
                               model_info,
                               model_fit_time, 
                               eigen_decomp_time = 0)
        
    }
    
    
    
    
    
    
}
