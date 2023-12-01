#!/bin/sh
#SBATCH --partition=short
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 0-0:10
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eberdan@hsph.harvard.edu

module load gcc/9.2.0 R/4.1.2
export R_LIBS_USER="~/R/4.1.2/library"

#Rscript diffbind_alpha_vs_beta.R 
#Rscript diffbind_alpha_vs_untreated.R
#Rscript diffbind_beta_vs_untreated.R
#Rscript diffbind_all.R
Rscript diffbind_all_dbpeaks.R 
