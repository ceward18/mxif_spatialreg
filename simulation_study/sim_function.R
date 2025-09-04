

# code to get tumor regions modified from scSpatialSIM
gen_one_image <- function(subjectID, imageID, sigma_spat, rho) {
    
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
    
    # add subject and image ID
    sim_dat$subjectID <- subjectID
    sim_dat$imageID <- imageID
    sim_dat$tumorID <- 1:nrow(sim_dat)
    
    # simulate spatial information using squared exponential
    sim_dat$e <- sim_spatial_mvn(sim_dat, sigma_spat, rho)
    
    sim_dat
}


# simulation function
sim_data_fn <- function(betas, nSub, nImagePerSub, 
                        rho, sigma_spat,
                        sigma_sub, sigma_image) {
    
    # subject specific random effect
    b_dat <- cbind.data.frame(subjectID = 1:nSub,
                              b = rnorm(nSub, sd = sigma_sub))
    
    # if only one image per subject, then don't need subject random effect
    if (nImagePerSub == 1) {
        b_dat$b <- 0
    }
    
    # image specific random effect
    c_dat <- cbind.data.frame(subjectID = rep(1:nSub, each = nImagePerSub),
                              imageID = 1:(nImagePerSub * nSub),
                              c = rnorm(nImagePerSub * nSub, sd = sigma_image))
    
    ### loop over images and generate data
    
    # first create data
    # x, y coordinates between 0 and 1 for each image
    imageList <- list()
    
    for (i in 1:nSub) {
        
        for (j in 1:nImagePerSub) {
            
            # sometimes doesn't work due to numerical issues
            dataOkay <- FALSE
            while(!dataOkay) {
                imageDat <- tryCatch(gen_one_image(subjectID = i,
                                                   imageID = (i-1)*nImagePerSub + j,
                                                   sigma_spat = sigma_spat, 
                                                   rho = rho),
                                     error = function(e) 'bad data')
                if (is.data.frame(imageDat)) {
                    dataOkay <- TRUE
                }
                
            }
            
            imageDat <- merge(imageDat, b_dat, by = 'subjectID', all.x = T)
            imageDat <- merge(imageDat, c_dat, by = c('subjectID', 'imageID'), all.x = T)
            
            # store data
            imageList[[(i-1)*nImagePerSub + j]] <- imageDat
        }
    }
    
    tumorDat <- do.call("rbind.data.frame", imageList)
    
    # half responder vs non responder
    nResponder <- ceiling(nSub/2)
    responderStatus <- data.frame(subjectID = 1:nSub,
                                  response = c(rep(1, nResponder),
                                               rep(0, nSub - nResponder)))
    tumorDat <- merge(tumorDat, responderStatus, by = 'subjectID', all.x = T)
    
    
    # probabilities 
    tumorDat$prob <- with(tumorDat, 
                          ilogit(betas[1] + betas[2] * response + b + c + e))
    
    # distribution of number of neighbors
    tumorDat$numNeighbors <- rnbinom(nrow(tumorDat), mu = 17, size = 15)
    
    # get outcome from binomial draw
    tumorDat$outcome <- rbinom(nrow(tumorDat), 
                               tumorDat$numNeighbors, 
                               tumorDat$prob)
    
    
    # select columns to output
    tumorDat[order(tumorDat$subjectID, tumorDat$imageID),
             c('subjectID', 'imageID', 'tumorID',
               'x', 'y', 'response', 'numNeighbors', 'outcome')]
    
}
