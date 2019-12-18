#!/bin/bash

#set the job name
#SBATCH --job-name=PrimateOrtho
#SBATCH -n 4
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03-00:00

#run the application
#PATHS

module load Python/3.6.6-foss-2018b
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/
mkdir ${OUTPUT}Primates


### Python script to create orthos
python3 PrimateQualOrthos.py ${OUTPUT}Primates/GeneLengthsQual.tsv
