#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=PC_sim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH --array=1-990
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
source /etc/profile
module load udunits/2.2.27.6_gcc11.3.0-rocky8
module load proj/7.2.1-rocky8
module load geos/3.7.1-rocky8
module load gdal/2.3.2-rocky8
module load openssl/3.1.1
module load R/4.4.0-openblas-rocky8

####### Run your script #########################
Rscript run_models.R $SLURM_ARRAY_TASK_ID
