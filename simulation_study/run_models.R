################################################################################
# Script to generate data, fit all possible models, and output useful summaries
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

library(parallel)
library(foreach)
library(MASS)
library(nimble)
library(spatstat.geom)
library(spatstat.random)
library(mvnfast)
library(Matrix)
library(spaMM)
library(doSNOW)

source('helper_functions.R')
source('sim_function.R')
source('fit_model.R')

# 396000
all_sims <- expand.grid(n_subjects = c(10, 30, 50),
                        n_image_sub = c(1, 5),
                        zero_distance = c(50, 100),
                        sigma_spat = c(0.5, 1, 2),
                        beta_idx = 1:11,
                        sim_number = 1:1000,
                        stringsAsFactors = F)

# prop outcome, beta1 corresponds to OR
ors <- c(1, 1.1, 1.25, 1.5, 2, 3)
ors <- c(rev(1/ors[-1]), ors)

beta_vals <- data.frame(beta_idx = 1:11,
                        beta_val = log(ors))

all_sims <- merge(all_sims, beta_vals, by = 'beta_idx', all.x = T)
all_sims <- all_sims[order(all_sims$n_subjects,
                           all_sims$n_image_sub,
                           all_sims$zero_distance,
                           all_sims$beta_idx,
                           all_sims$sim_number),]
rownames(all_sims) <- NULL


# do batches
batch_size <- 400
batch_idx <- batch_size * (idx - 1) + 1:batch_size

for (i in batch_idx) {
    
    # get specifications
    n_subjects <- all_sims$n_subjects[i]
    n_image_sub <- all_sims$n_image_sub[i]
    zero_distance <- all_sims$zero_distance[i]
    sigma_spat <- all_sims$sigma_spat[i]
    beta_val <- all_sims$beta_val[i]
    sim_number <- all_sims$sim_number[i]
    beta_idx <- all_sims$beta_idx[i]
    
    print(paste0('n sub: ', n_subjects, 
                 ', n images per sub: ', n_image_sub, 
                 ', zero distance: ', zero_distance, 
                 ', sigma^2: ', sigma_spat, 
                 ', beta: ', round(beta_val, 2),
                 ', simulation: ', sim_number))
    
    ############################################################################
    
    ### simulate data
    
    # intercept specifies proportion for non-responders
    beta0 <- logit(0.05)
    
    # rho depends on zero_distance
    rho <- optim(0, function(x) (exp(-zero_distance^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par
    
    set.seed(i)
    
    # take up to 4-5 min for 500 image datasets (5 sec for small data)
    sim_data <- sim_data_fn(betas = c(beta0, beta_val), 
                            nSub = n_subjects, 
                            nImagePerSub = n_image_sub, 
                            rho = rho,
                            sigma_spat = sigma_spat,    # spatial SD
                            sigma_sub = 0.4,            # between subjects SD
                            sigma_image = 0.4)          # between images SD
    
    ############################################################################
    # fit models
    
    model_info <- data.frame(n_subjects, n_image_sub, 
                             zero_distance, sigma_spat, sim_number, beta_idx) 
    
    # fit no_corr and pc models in parallel
    model_types <- c('no_corr', 
                     'pc_exp', 
                     'pc_sqexp')
    
    res_models <- list()
    
    # takes about 11 min on largest datasets 
    
    for (j in 1:length(model_types)) {
        
        print(model_types[j])
        
        res_models[[j]] <- fit_model(model_type = model_types[j],
                                     sim_data, model_info, beta0, beta_val)
    }
    
    all_results <- do.call("rbind.data.frame", res_models)
    all_results$n_cells <- nrow(sim_data)
    rownames(all_results)  <- NULL
    
    # save batch results
    if (i == batch_idx[1]) {
        all_batches <- all_results
    } else {
        all_batches <- rbind.data.frame(all_batches, all_results)
    }
    
    
}

# save output in RDS form
saveRDS(all_batches, paste0('./output/res_batch_', sprintf("%03d",idx), '.rds'))


