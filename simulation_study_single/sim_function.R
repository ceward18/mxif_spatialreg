

# code to get tumor regions modified from scSpatialSIM
gen_one_image <- function(sigma_spat, rho) {
    
    # simulate a point process with lambda == 0.008
    get_pp <- rpoispp(0.008, win = as.owin(c(0, 750,
                                             0, 750)))
    sim_dat <- data.frame(x = get_pp$x, y = get_pp$y)
    
    #simulate centers centers of groups
    k <- 5
    gauss_tab <- data.frame(x = runif(k, 0, 750), 
                            y = runif(k, 0, 750),
                            sd.x = runif(k, 5, 300), 
                            sd.y = runif(k, 5, 300))
    
    #spearman correlations for how the distribution falls in 2D space
    gauss_tab$rho <- runif(k,
                           -gauss_tab$sd.x * gauss_tab$sd.y, 
                           gauss_tab$sd.x * gauss_tab$sd.y)
    
    # compute probability from each kernel and then take max
    prob_kernel <- vapply(1:k, function(a){
        diff <- cbind((sim_dat$x - gauss_tab$x[a]), (sim_dat$y - gauss_tab$y[a]))
        sigma_mat <- matrix(c(gauss_tab$sd.x[a]^2, rep(gauss_tab$rho[a], 2),
                              gauss_tab$sd.y[a]^2), nrow = 2, ncol = 2)
        z <- diag(diff %*% solve(sigma_mat) %*% t(diff))
        exp(-z) 
    }, FUN.VALUE = numeric(nrow(sim_dat)))
    sim_dat$prob <- apply(prob_kernel, 1, max) * 0.9
    sim_dat$tumor_cell <- rbinom(nrow(sim_dat), 1, sim_dat$prob)
    sim_dat <- sim_dat[which(sim_dat$tumor_cell == 1), c('x', 'y')]
    
    sim_dat$tumorID <- 1:nrow(sim_dat)
    
    # simulate spatial information using squared exponential
    sim_dat$e <- sim_spatial_mvn(sim_dat, sigma_spat, rho)
    
    sim_dat
}


# simulation function for single image with binary cell-level covariate
sim_data_fn <- function(betas, rho, sigma_spat) {
    
    ### loop over images and generate data
    
    # sometimes doesn't work due to numerical issues
    dataOkay <- FALSE
    while(!dataOkay) {
        tumorDat <- tryCatch(gen_one_image(sigma_spat = sigma_spat, 
                                           rho = rho),
                             error = function(e) 'bad data')
        if (is.data.frame(tumorDat)) {
            dataOkay <- TRUE
        }
        
    }
    
    # ~50% of cells have binary marker
    tumorDat$marker <- rbinom(nrow(tumorDat), 1, 0.5)
    
    # probabilities 
    tumorDat$prob <- with(tumorDat, 
                          ilogit(betas[1] + betas[2] * marker + e))
    
    # distribution of number of neighbors
    tumorDat$numNeighbors <- rpois(nrow(tumorDat), 5)
    
    # get outcome from binomial draw
    tumorDat$outcome <- rbinom(nrow(tumorDat), 
                               tumorDat$numNeighbors, 
                               tumorDat$prob)
    
    
    # select columns to output
    tumorDat[order(tumorDat$tumorID),
             c('tumorID', 'x', 'y', 'marker', 'numNeighbors', 'outcome')]
    
}
