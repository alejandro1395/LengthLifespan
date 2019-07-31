#!/bin/bash

#set the job name
#SBATCH --job-name=PrimateOrtho
#SBATCH -n 4
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00

#run the application
#PATHS

OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/
mkdir ${OUTPUT}Primates

module load Python


### Python script to create orthos
python3 PrimateOrthoKeep.py ${OUTPUT}Primates/ortho1to1.tsv
