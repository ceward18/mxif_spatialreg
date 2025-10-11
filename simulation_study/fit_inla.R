################################################################################
# INLA implementation


fit_inla <- function(sim_data, n_image_sub,  beta_val, beta2, model_info) {
  
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
  
  
  mesh_image <- inla.mesh.2d(
    loc.domain = cbind(
      x = c(0, 750, 750, 0),
      y = c(0, 0, 750, 750)
    ), 
    max.edge = c(30, 100),
    offset = c(50, 150))
  # plot(mesh_image)
  
  
  # SPDE
  spde <- inla.spde2.pcmatern(mesh = mesh_image, 
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
  
  # half normal prior for sds
  HN_prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);  
"
  
  
  # fit model
  
  if (n_image_sub == 1) {
    
    formulae <- y ~ 0 + b0 + response + fov_type + 
      f(imageID, model = "iid", hyper = list(prec = list(prior = HN_prior, initial = 1))) +
      f(i, model = spde,
        group = i.group, control.group = list(model = 'iid')) 
    
  } else {
    
    formulae <- y ~ 0 + b0 + response + fov_type + 
      f(imageID, model = "iid", hyper = list(prec = list(prior = HN_prior, initial = 1))) +
      f(subjectID, model = "iid", hyper = list(prec = list(prior = HN_prior, initial = 1))) +
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
              control.compute=list(config = TRUE),
              num.threads = 10)
  
  end_model <- Sys.time()
  
  model_fit_time <-  as.numeric(end_model - start_model, units = 'secs')
  
  get_summary_table_inla(res, beta_val, beta2, 
                         model_info, model_fit_time)
  
  
}



