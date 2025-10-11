#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=PC_sim_inla
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --array=1-17
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
source /etc/profile
module load udunits/2.2.27.6_gcc11.3.0-rocky8
module load proj/7.2.1-rocky8
module load geos/3.7.1-rocky8
module load gdal/2.3.2-rocky8
module load openssl/3.1.1
module load gcc/13.1.0-5z64cho
module load R/4.4.0-openblas-rocky8

####### Run your script #########################
Rscript run_inla_model.R $SLURM_ARRAY_TASK_ID
