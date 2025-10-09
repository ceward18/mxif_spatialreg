################################################################################
# INLA implementation


fit_inla <- function(sim_data, n_image_sub,  beta_val, beta2) {
    
    # remove those with no neighbors
    sub_dat <- sim_data[sim_data$numNeighbors > 0,]
    
    # total number of images (for grouping)
    k <- length(unique(sub_dat$imageID))
    
    
    start_model <- Sys.time()
    
    ### create mesh
    # max edge controls triangles (lower = more triangles)
    # offset is inner and outer ring
    # cutoff is ?
    
    # inner triangle edge ~ range / 3 to range / 5 (for range = 100, 20-33)
    # outer 2-4x the inner max edge ()
    
    
    mesh_image <- fmesher::fm_mesh_2d_inla(
        loc.domain = cbind(
            x = c(0, 750, 750, 0),
            y = c(0, 0, 750, 750)
        ), 
        max.edge = c(30, 100),
        offset = c(50, 150))
    plot(mesh_image)
    
    
    # SPDE
    spde <- INLA::inla.spde2.pcmatern(mesh = mesh_image, 
                                      prior.range = c(200, 0.01), # P(range < 200) = 0.01
                                      prior.sigma = c(0.2, 0.01)) # P(sigma > 0.2) = 0.01
    
    # group indexing
    iset <- INLA::inla.spde.make.index('i', 
                                       n.spde = spde$n.spde,
                                       n.group = k)
    
    # projection matrix
    A <- INLA::inla.spde.make.A(mesh = mesh_image,
                                loc = cbind(sub_dat$x, sub_dat$y), 
                                group = sub_dat$imageID) 
    
    # stack
    sdat <- INLA::inla.stack(
        data = list(y = sub_dat$outcome, n = sub_dat$numNeighbors), 
        A = list(A, 1), 
        effects = list(iset,
                       data.frame(b0 = rep(1, nrow(sub_dat)),
                                  response = sub_dat$response, 
                                  fov_type = sub_dat$fov_type,
                                  subjectID = sub_dat$subjectID,
                                  imageID = sub_dat$imageID)
        ),
        tag = 'stdata') 
    
    
    # fit model
    
    if (n_image_sub == 1) {
        
        formulae <- y ~ 0 + b0 + response + fov_type + 
            f(imageID, model = "iid") +
            f(i, model = spde,
              group = i.group, control.group = list(model = 'iid')) 
        
    } else {
        
        formulae <- y ~ 0 + b0 + response + fov_type + 
            f(imageID, model = "iid") +
            f(subjectID, model = "iid") +
            f(i, model = spde,
              group = i.group, control.group = list(model = 'iid')) 
    }
    
    
    
    
    
    
    # Model fitting
    res <- inla(formulae, 
                data = inla.stack.data(sdat), 
                family = "binomial",
                Ntrials = n,
                control.predictor = list(compute = TRUE,
                                         A = inla.stack.A(sdat),
                                         link = 1),
                control.compute=list(config = TRUE))
    
    end_model <- Sys.time()
    
    model_fit_time <-  as.numeric(end_model - start_model, units = 'secs')
    
    get_summary_table_inla(res, beta_val, beta2, sigma_spat_use, 
                           model_info, model_fit_time)
    
    
}




res$summary.fixed
res$summary.hyperpar

res$summary.random$imageID
res$summary.random$subjectID

length(res$marginals.hyper)
par(mfrow = c(2, 2))
for (j in 1:length(res$marginals.hyper)) {
    plot(res$marginals.hyper[[j]], type = 'l', 
         xlab = names(res$marginals.hyper)[j], ylab = 'Density')
    abline(v = c(1 /(0.4^2), 1 /(0.4^2), 100, 2)[j],
           col = 2)
}



samples <- inla.posterior.sample(5000, res)


# Example: posterior of exp(beta_response)
marg_beta_response <- res$marginals.fixed$"response"

# Transform: exp(beta)
marg_OR_response <- inla.tmarginal(function(x) exp(x), marg_beta_response)

# Summarize the posterior for OR
inla.zmarginal(marg_OR_response)


plot(marg_OR_response, type = "l",
     main = "Posterior of exp(beta_response)",
     xlab = "Odds Ratio", ylab = "Density")





