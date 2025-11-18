#!/bin/bash

#SBATCH --job-name 05_plotting
#SBATCH -n 10                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 05_plot_looping_regions_for_visualScreening.out # File to which STDOUT will be written
#SBATCH -e 05_plot_looping_regions_for_visualScreening.out # File to which STDERR will be written

conditions="larvae_DWT"
assay=hic

for condition in $conditions ;
#for condition in $1 ;
do
    for file in $(ls -1 Loops*${assay}*${condition}*5_P*FDR_resD*merged.bedpe);
    do
	echo $file
	
	mkdir -p ${file%.bedpe}

	sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n $file > _tmp_${file} ; mv _tmp_${file} ${file}
	
	python /work/user/mdistefano/2022_08_22_microC_eyeDiscs_Bernd/04_loopCalling_analysis/loopDetection_analysis/mustache_loopcaller/scripts_clean/05_plot_looping_regions_for_visualScreening.py ${file%.bedpe} ${condition}
	#exit
    done
done
	    
