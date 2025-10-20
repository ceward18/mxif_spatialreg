# mxif_spatialreg

This repository contains code associated with the simulation study and data analysis presented in "Scalable Inference for Spatial Generalized Mixed Models for use in Multiplexed Imaging Analysis" by Ward et al. (2025).

The provided code implements the proposed eigen-decomposition approach for fitting the spatial GLMM for co-localization, as well as the SPDE approximation using R-INLA. 


## simulation_study

/output/ - this directory contains batches of stored results from model fitting on the HPC system.

/results/ - this directory contains processed results concatenated across all batches.

combine_batches.R - this R script takes the batched results from the /output/ directory, concatenates and processes them, and then outputs the processed files to the /results/ directory.

create_figures.R - this R script recreates Figures 1-4 in the manuscript.

fit_inla.R - this R script contain one function that fits the SPDE model with R-INLA to a single simulated dataset and returns relevant output.

fit_model.R - this R script contains one function that fits the eigen-decomposition and no spatial correlation models to a single simulated dataset and returns relevant output.

helper_functions.R - this R script contains helper functions for simulating data, model fitting, and summarizing model output.

run_inla_model.R - this R script runs the subset of simulated datasets where the SPDE approach is fitted using R-INLA.

run_models.R - this R script runs the full simulations for the eigen-decomposition approach and no spatial correlation models.

sim_functions.R - this R script contains functions to generate spatial locations on index cells and simulate the co-localization outcome for the given study design.




