#!/bin/sh
#SBATCH --partition=priority
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 2-0:00
#SBATCH -c 20
#SBATCH --mem=70G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eberdan@hsph.harvard.edu


bcbio_nextgen.py ../config/mur_mel.yaml -n 20
