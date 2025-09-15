################################################################################
# helper functions
################################################################################

sim_spatial_mvn <- function(sim_dat, sigma_spat, rho) {
    # get distance matrix
    coords <- cbind(sim_dat$x, sim_dat$y)
    dist_mat <- as.matrix(dist(coords))
    # compute correlation
    cov_mat <- sigma_spat^2 * exp(- dist_mat^2 /(2 * rho^2))
    # use mvnfast package
    c(rmvn(1, rep(0, nrow(sim_dat)), cov_mat))
}


# get correlation matrix for a single image 
get_correlation <- function(fov_data, zero_distance, sigma_spat, cor_type, induce_sparsity) {
    # get distance matrix
    coords <- cbind(fov_data$x, fov_data$y)
    dist_mat <- as.matrix(dist(coords))
    
    # convert to spatial correlation
    # find rho value for which correlation at zero_distance = 0.001
    if (cor_type == 'exp') {
        # corr = exp(-d / rho)
        
        rho <- optim(0, function(x) (exp(-zero_distance / x) - 0.001)^2, 
                     method = 'Brent',
                     lower = 0, upper = 500)$par
        cor_mat <- sigma_spat^2 * exp(-dist_mat / rho)
        
    } else if (cor_type == 'sqexp') {
        # corr = exp(-d^2 / (2 * rho^2))
        
        rho <- optim(0, function(x) (exp(-zero_distance^2 / (2 * x^2)) - 0.001)^2, 
                     method = 'Brent',
                     lower = 0, upper = 500)$par
        
        cor_mat <- sigma_spat^2 * exp(- dist_mat^2 /(2 * rho^2))
    }
    
    if (induce_sparsity) {
        cor_mat[cor_mat < 0.001] <- 0
    }
    
    cor_mat
}

# get estimated zero_distance and sigma_spat using variogram
variogram_est <- function(fov_data) {
    
    # variance is on the scale of logit(p)
    fov_data <- fov_data[fov_data$numNeighbors > 0,]
    fov_data$prop <- fov_data$outcome / fov_data$numNeighbors
    fov_data$logit_prop <- log(fov_data$prop / (1 - fov_data$prop)) 
    fov_data$logit_prop[fov_data$prop == 0] <- log(0.01 / 0.99) 
    fov_data$logit_prop[fov_data$prop == 1] <- log(0.99 / 0.01) 
    coordinates(fov_data) = ~x+y
    var_out <- variogram(logit_prop~1, fov_data, width = 2)
    
    
    fit_model <- fit.variogram(var_out, model=vgm('Sph'))
    
    c('sigma_spat' = fit_model$psill[2],
      'zero_distance' = fit_model$range[2])
    
} 

get_summary_table <- function(spfit, sim_data, 
                              beta_val, beta2, 
                              sigma_spat_use, model_type, model_info,
                              model_fit_time, eigen_decomp_time) {
    
    
    start <- Sys.time()
    
    ci <- confint(spfit)
    
    # exponentiate to OR scale
    
    sum_tab <- data.frame(coef = c('OR_R', 'OR_FOV'),
                          est = exp(ci[2:3,3]),
                          se = sqrt(diag(vcov(spfit)$cond)[2:3]),
                          lower = exp(ci[2:3,1]),
                          upper = exp(ci[2:3,2]),
                          truth = c(exp(beta_val), beta2),
                          model_info,
                          model_type = model_type)
    
    
    end <- Sys.time()
    
    summary_time <- as.numeric(end - start, units = 'secs')
    
    
    if (model_type == 'no_corr') {
        sum_tab$sigma_spat_est <- NA
    } else {
        sum_tab$sigma_spat_est <- sigma_spat_use
    }
    
    sum_tab$eigen_decomp_time <- eigen_decomp_time
    sum_tab$model_fit_time <- model_fit_time 
    sum_tab$summary_time <- summary_time
    sum_tab$time <- eigen_decomp_time + model_fit_time + summary_time
    
    
    sum_tab
}

