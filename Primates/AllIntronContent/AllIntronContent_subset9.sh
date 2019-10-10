#!/bin/bash

set the job name
#SBATCH --job-name=AllIntronContent
#SBATCH -n 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03-00:00

#run the application
#PATHS

module load foss
module load Python

	python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/AllIntronContent.py /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/ortho1to1_subset9.tsv /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/AllIntronContent_subset9.tsv > /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/AllIntronContent/out/AllIntronContent_subset9.out
