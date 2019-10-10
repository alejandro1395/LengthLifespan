#!/bin/bash

set the job name
#SBATCH --job-name=CDSfromGene
#SBATCH -n 1
#SBATCH -o /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene/out/CDSfromGene_subset704.out
#SBATCH -e /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene/qu/CDSfromGene_subset704.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03-00:00

#run the application
#PATHS

module load foss
module load Python

        python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene.py /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/ortho1to1_subset704.tsv /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/CDSfromGene_subset704.tsv
