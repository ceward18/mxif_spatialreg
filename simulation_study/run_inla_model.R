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
library(fmesher)

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
                        sim_number = 1:20,
                        stringsAsFactors = F)


all_sims <- merge(all_sims, beta_vals, by = 'beta_idx', all.x = T)

all_sims <- all_sims[order(all_sims$n_subjects,
                           all_sims$n_image_sub,
                           all_sims$beta_idx,
                           all_sims$sim_number),]
rownames(all_sims) <- NULL

tmp <- all_sims[seq(1, nrow(all_sims), 5),]
rownames(tmp) <- NULL

# do batches
batch_size <- 200
batch_idx <- batch_size * (idx - 1) + 1:batch_size


for (i in nrow(all_sims)) {
    
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
    
    all_results <- fit_inla(sim_data, n_image_sub,  beta_val, beta2)
    
    
    
    
    
    
    # fit the model using exponential
    start_model <- Sys.time()
    
    if (n_image_sub == 1) {
        # no subject random effect
        spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ response + fov_type +
                           (1|imageID) +
                           MaternIMRFa(1|x+y %in% imageID),
                       data=sub_dat,
                       family = binomial,
                       fixed=list(nu=1),
                       method = 'REML')
    } else {
        spfit <- fitme(cbind(outcome, numNeighbors-outcome) ~ response + fov_type +
                           (1|subjectID) +
                           (1|imageID) +
                           Matern(1|x+y %in% imageID),
                       data=sub_dat,
                       family = binomial,
                       fixed=list(nu=100),
                       method = 'REML')
    }

    end_model <- Sys.time()
    
    
    
    
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
saveRDS(all_batches, paste0('./output/res_batch_', sprintf("%04d",idx), '.rds'))


