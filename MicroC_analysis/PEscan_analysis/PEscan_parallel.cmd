#!/bin/bash
#SBATCH --job-name peScan
#SBATCH -n 12                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=16G               # allocated memory
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00             # Runtime in D-HH:MM
#SBATCH -o /work/cavalli/2022_06_08_Project_on_PREs_contacts/peScan_analysis/PEscan.out   # File to which STDOUT will be written
#SBATCH -e /work/cavalli/2022_06_08_Project_on_PREs_contacts/peScan_analysis/PEscan.out   # File to which STDERR will be written

#for sample in PH18 PH29 ;
#for sample in Embryo WD LD PH18 ;
for sample in Embryo ;
do
    Rscript ./scripts/PEscan_parallel.R $sample &>> PEscan_parallel_${sample}.out &
done
