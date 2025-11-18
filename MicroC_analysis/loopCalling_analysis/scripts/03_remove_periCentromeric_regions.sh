#!/bin/bash

#SBATCH --job-name 03_step
#SBATCH -n 10                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 03_remove_periCentromeric_regions.out # File to which STDOUT will be written
#SBATCH -e 03_remove_periCentromeric_regions.out # File to which STDERR will be written

assembly=dm6
blackListFile=${PWD}/scripts/${assembly}_regions_to_exclude.bed
conditions="larvae_DWT"
resolutions="1000 2000 4000 8000 10000 20000 40000"

for dir in $(ls -1 | grep calledLoops_ | grep calledLoops_8000_upTo_1348131bp) ;	   
do
    cd $dir
    echo $dir

    for file in $(ls -1 hic*_sparsityThreshold_1.00_*_5) ;
    do
	outFile=${file}_PeriCentromere_filtered
	if [[ ! -e ${outFile} ]];
	then
	    touch ${outFile}
	    echo "${file} -> ${outFile}"
	    echo "Initial number of loops in $(wc -l ${file} | awk '{print $1-1}')"	    

	    echo "Filtering-out loops in peri-centromeric regions"
	    
	    echo "Remove loops in blacklisted regions (e.g., pericentromeric regions)"
	    awk '{if($4=="BL"){n++;chrom[n]=$1;start[n]=$2;end[n]=$3;next}else{flag=0;for(i=1;i<=n;i++){if($1==chrom[i]){for(j=start[i];j<=end[i];j+=100){if($2<=j && j<=$3){flag=1}}};if(flag==0 && $4==chrom[i]){for(j=start[i];j<=end[i];j+=10){if($5<=j && j<=$6){flag=1; next}}}}; if(flag==0){print $0}}}' <(awk '{print $1,$2,$3,"BL"}' ${blackListFile}) ${file} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | grep -v BIN > ${outFile}

	    echo "Number of loops after filtering blacklisted regions $(wc -l ${outFile} | awk '{print $1-1}')"
	    echo ""
	fi
    done # Close cycle over $file
    cd .. # Exit $dir
    
done # Close cycle over $dir
