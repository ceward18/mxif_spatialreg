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
get_correlation <- function(fov_data, zero_distance, cor_type) {
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
        cor_mat <- exp(-dist_mat / rho)
        
    } else if (cor_type == 'sqexp') {
        # corr = exp(-d^2 / (2 * rho^2))
        
        rho <- optim(0, function(x) (exp(-zero_distance^2 / (2 * x^2)) - 0.001)^2, 
                     method = 'Brent',
                     lower = 0, upper = 500)$par
        
        cor_mat <- exp(- dist_mat^2 /(2 * rho^2))
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
    
    c('sigma_spat' = sqrt(fit_model$psill[2]),
      'zero_distance' = fit_model$range[2])
    
} 

# general function to get confidence interval for beta1
get_confint <- function(spfit) {
    
    betaCovMat <- vcov(spfit)
    betaHat <- spfit$fixef
    
    # response
    LVec <- c(0,1)
    seLogitScale <- as.vector(sqrt(t(LVec) %*% betaCovMat %*% LVec))
    ci2 <- c(LVec %*% betaHat)  + c(-1.96, 1.96) * seLogitScale
    
    matrix(exp(ci2), nrow = 1, byrow = T)
    
}


get_summary_table <- function(spfit, sim_data, 
                              beta0, beta_val,
                              sigma_spat_use, zero_dist_use,
                              model_type, model_info,
                              model_fit_time, eigen_decomp_time) {
    
    
    start <- Sys.time()
    
    ci <- get_confint(spfit)
    
    # exponentiate to OR scale
    # also want probabilities of outcome for responder and non-responders
    prob_dat <- data.frame(marker = c(0, 1))
    if (grepl('pc_', model_type)) {
        prob_dat$spatial_offset = mean(unique(sim_data$spatial_offset))
    }
    
    prob_dat <- cbind.data.frame(prob_dat,
                                 est = predict(spfit, 
                                               newdata = prob_dat, 
                                               re.form = NA,
                                               type = 'response'))
    # compute variance on scale of linear predictor (not on probability scale!)
    prob_dat$var <- get_fixefVar(spfit,
                                 newdata=prob_dat, 
                                 re.form = NA)
    prob_dat$lower <- ilogit(logit(prob_dat$est) - 1.96 * sqrt(prob_dat$var))
    prob_dat$upper <- ilogit(logit(prob_dat$est) + 1.96 * sqrt(prob_dat$var))
    
    sum_tab <- data.frame(coef = c('OR', 'p_0', 'p_1'),
                          est = c(exp(spfit$fixef)[2], prob_dat$est),
                          se = c(sqrt(diag(vcov(spfit)))[2], sqrt(prob_dat$var)),
                          lower = c(ci[1,1], prob_dat$lower),
                          upper = c(ci[1,2], prob_dat$upper),
                          truth = c(exp(beta_val), ilogit(c(beta0, beta0 + beta_val))),
                          model_info,
                          model_type = model_type,
                          sigma_spat_est = sigma_spat_use,
                          zero_dist_est = zero_dist_use)
    
    
    end <- Sys.time()
    
    summary_time <- as.numeric(end - start, units = 'secs')

    
    sum_tab$eigen_decomp_time <- eigen_decomp_time
    sum_tab$model_fit_time <- model_fit_time 
    sum_tab$summary_time <- summary_time
    sum_tab$time <- eigen_decomp_time + model_fit_time + summary_time
    
    
    sum_tab
}



get_summary_table_nngp <- function(nngp_mod, 
                                   beta0, beta_val,
                                   model_type, model_info,
                                   model_fit_time, eigen_decomp_time) {
    
    
    start <- Sys.time()
    
    # need 95% CI for OR and probabilities
    p0 <- ilogit(nngp_mod$p.beta.samples[nngp_mod$sub.sample$start:nngp_mod$sub.sample$end, 1])
    p1 <- ilogit(rowSums(nngp_mod$p.beta.samples[nngp_mod$sub.sample$start:nngp_mod$sub.sample$end, ]))
    
    or <- (p1 / (1 - p1)) / (p0 / (1 - p0))
    
    sum_tab <- data.frame(coef = c('OR', 'p_0', 'p_1'),
                          est = c(mean(or), mean(p0), mean(p1)),
                          se = c(sd(or), sd(p0), sd(p1)),
                          lower = c(quantile(or, probs = 0.025), 
                                    quantile(p0, probs = 0.025), 
                                    quantile(p1, probs = 0.025)),
                          upper = c(quantile(or, probs = 0.975), 
                                    quantile(p0, probs = 0.975), 
                                    quantile(p1, probs = 0.975)),
                          truth = c(exp(beta_val), ilogit(c(beta0, beta0 + beta_val))),
                          model_info,
                          model_type = model_type,
                          sigma_spat_est = median(nngp_mod$p.theta.samples[nngp_mod$sub.sample$start:nngp_mod$sub.sample$end, 1]),
                          zero_dist_est = median(sqrt(-log(0.001)) / nngp_mod$p.theta.samples[nngp_mod$sub.sample$start:nngp_mod$sub.sample$end, 3]))
    
    
    end <- Sys.time()
    
    summary_time <- as.numeric(end - start, units = 'secs')
    
    sum_tab$eigen_decomp_time <- eigen_decomp_time
    sum_tab$model_fit_time <- model_fit_time 
    sum_tab$summary_time <- summary_time
    sum_tab$time <- eigen_decomp_time + model_fit_time + summary_time
    
    
    sum_tab
}





