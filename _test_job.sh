#!/bin/bash
#SBATCH --job-name=parallel_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shr264@ufl.edu
#SBATCH --account=statistics
#SBATCH --qos=statistics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8gb
#SBATCH --time=48:00:00
#SBATCH --output=_test_%j.out
date;hostname;pwd

module load intel/2016.0.109 R/3.3.0

Rscript _test.R

date
