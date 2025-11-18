#!/bin/bash

#SBATCH --job-name 07_recentering
##SBATCH -n 10                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 150Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 07_recentering.out # File to which STDOUT will be written
#SBATCH -e 07_recentering.out # File to which STDERR will be written

Rscript ../scripts_clean/07_recenter_loops_by_maximum_shamanScore.R $1
