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

OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/
SRC=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/
mkdir -p ${SRC}CDSfromGene

module load Python

for i in {1..721};
do FILE=${SRC}CDSfromGene/CDSfromGene_subset${i}.sh
cat  > $FILE << EOT
#!/bin/bash

set the job name
#SBATCH --job-name=CDSfromGene
#SBATCH -n 1
#SBATCH -o /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene/out/CDSfromGene_subset${i}.out
#SBATCH -e /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene/qu/CDSfromGene_subset${i}.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03-00:00

#run the application
#PATHS

module load foss
module load Python

        python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/src/CDSfromGene.py /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/ortho1to1_subset${i}.tsv /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/RegLifespanJuly2019_PhD/LengthGenes_July2019/results/Primates/CDSfromGene_subset${i}.tsv
EOT

sbatch $FILE

done
