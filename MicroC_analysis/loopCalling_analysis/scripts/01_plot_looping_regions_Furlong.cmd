#!/bin/bash

#SBATCH --job-name 01_plot
#SBATCH -n 10                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 150Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 01_plot_looping_regions_Furlong.out # File to which STDOUT will be written
#SBATCH -e 01_plot_looping_regions_Furlong.out # File to which STDERR will be written 

python scripts/01_plot_looping_regions_Furlong.py $1 $2
