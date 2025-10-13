################################################################################
# Script to generate data, fit all possible models, and output useful summaries
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

library(MASS)
library(nimble)
library(spatstat.geom)
library(spatstat.random)
library(mvnfast)
library(Matrix)
library(INLA)

source('helper_functions.R')
source('sim_function.R')
source('fit_inla.R')

# prop outcome, beta1 corresponds to OR
ors <- c(1, 1.25, 1.5, 2, 3)
ors <- c(rev(1/ors[-1]), ors)

beta_vals <- data.frame(beta_idx = 1:length(ors),
                        beta_val = log(ors))

# 3240
all_sims <- expand.grid(n_subjects = c(10, 30, 50),
                        n_image_sub = c(1, 5),
                        sigma_spat = c(0.5, 1, 2),
                        beta_idx = 1:nrow(beta_vals),
                        sim_number = 1:100,
                        stringsAsFactors = F)

# 16200
all_sims <- merge(all_sims, beta_vals, by = 'beta_idx', all.x = T)

all_sims <- all_sims[order(all_sims$n_image_sub,
                           all_sims$n_subjects,
                           all_sims$beta_idx,
                           all_sims$sim_number),]
rownames(all_sims) <- NULL


# do batches
# batch 1-60 get 1:8100 (all nimages = 1)
# batches 8101:10800 is 10/5 ~ 10 min
# batches 10801:13500 is 30/5 ~ 20 min (need high memory)
# batches 13501:16200 is 50/5 ~ 40 min (need high memory)

if (idx <= 60) {
    batch_size <- 135
    batch_idx <- batch_size * (idx - 1) + 1:batch_size
} else if (idx <= 85){ # idx 61-85, 25 batches of 108
    batch_size <- 108
    batch_idx <- (135 * 60) + 108 * (idx - 60 - 1) + 1:batch_size
} else if (idx <= 135){ # idx 86-135, 50 batches of 54
    batch_size <- 54
    batch_idx <- (135 * 60) + (108 * 25) + 54 * (idx - 85 - 1) + 1:batch_size
} else {               # idx 136-235, 100 batches of 27
    batch_size <- 27
    batch_idx <- (135 * 60) + (108 * 25) + (54 * 50) + 27 * (idx - 135 - 1) + 1:batch_size
}
batch_idx

for (i in batch_idx) {
    
    # get specifications
    n_subjects <- all_sims$n_subjects[i]
    n_image_sub <- all_sims$n_image_sub[i]
    sigma_spat <- all_sims$sigma_spat[i]
    beta_val <- all_sims$beta_val[i]
    beta_idx <- all_sims$beta_idx[i]
    sim_number <- all_sims$sim_number[i]
    
    print(paste0('n sub: ', n_subjects, 
                 ', n images per sub: ', n_image_sub, 
                 ', sigma^2: ', sigma_spat, 
                 ', beta: ', round(beta_val, 2),
                 ', simulation: ', sim_number))
    
    ############################################################################
    
    ### simulate data
    
    # intercept specifies proportion for non-responders in the core
    beta0 <- logit(0.05)
    
    # beta2 is difference between core and interface, fixed at OR 1.5
    beta2 <- log(1.5)
    
    # rho depends on zero_distance, fixed at 100
    rho <- optim(0, function(x) (exp(-100^2 / (2 * x^2)) - 0.001)^2, 
                 method = 'Brent',
                 lower = 0, upper = 500)$par
    
    set.seed(i)
    
    # take up to 4-5 min for 500 image datasets (5 sec for small data)
    sim_data <- sim_data_fn(betas = c(beta0, beta_val, beta2), 
                            nSub = n_subjects, 
                            nImagePerSub = n_image_sub, 
                            rho = rho,
                            sigma_spat = sigma_spat,    # spatial SD
                            sigma_sub = 0.4,            # between subjects SD
                            sigma_image = 0.4)          # between images SD
    
    ############################################################################
    # fit models
    
    model_info <- data.frame(n_subjects, n_image_sub, sigma_spat, sim_number, beta_idx) 
    
    
    # fit model and extract estimates + 95% CIs
    # 10, 1 = 41 sec
    # 10, 5 = 377 sec (6 min)
    # 30, 1 = 149 sec (2 min)
    # 30, 5 = 1084 sec (18 min)
    # 50, 1 = 192.1932 (3 min)
    # 50, 5 = 
    
    all_results <- fit_inla(sim_data, n_image_sub,  beta_val, beta2, model_info)
    
    
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
saveRDS(all_batches, paste0('./output/res_inla_batch_', sprintf("%03d",idx), '.rds'))


