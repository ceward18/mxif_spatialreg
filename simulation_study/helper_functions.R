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
    
    c('sigma_spat' = sqrt(fit_model$psill[2]),
      'zero_distance' = fit_model$range[2])
    
} 

get_confint <- function(spfit) {
    
    betaCovMat <- vcov(spfit)
    betaHat <- spfit$fixef
    
    # response
    LVec <- c(0,1)
    seLogitScale <- as.vector(sqrt(t(LVec) %*% betaCovMat %*% LVec))
    ci2 <- c(LVec %*% betaHat)  + c(-1.96, 1.96) * seLogitScale
    
    matrix(exp(ci2), nrow = 1, byrow = T)
    
}


get_confint_1 <- function(spfit) {
    
    betaCovMat <- vcov(spfit)
    betaHat <- spfit$fixef
    
    # estimated proportion
    
    LVec <- 1
    seLogitScale <- as.vector(sqrt(t(LVec) %*% betaCovMat %*% LVec))
    ci <- c(LVec %*% betaHat)  + c(-1.96, 1.96) * seLogitScale
    
    matrix(exp(ci), nrow = 1, byrow = T)
    
    
}

get_summary_table <- function(spfit, sim_data, 
                              beta0, beta_val,
                              sigma_spat_use, model_type, model_info,
                              model_fit_time, eigen_decomp_time) {
    
    
    start <- Sys.time()
    
    ci <- get_confint(spfit)
    
    # exponentiate to OR scale
    # also want probabilities of outcome for responder and non-responders
    prob_dat <- data.frame(response = c(0, 1))
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

