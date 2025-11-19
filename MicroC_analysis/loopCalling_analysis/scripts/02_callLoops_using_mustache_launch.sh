#!/bin/bash

conditions="Embryo_WT WD_WT ED_PH18"
resolutions="1000 2000" ; maxDists=400000
#resolutions="1000 2000 4000" ; maxDists=800000
#resolutions="4000 8000 10000 20000" ; maxDists=3200000
#resolutions="20000 40000" ; maxDists=33000000
sparsityThresholds="1.00"
sigmaZeros=$(seq 3.6 -0.1 0.6)
chroms="chr2L chr2R chr3L chr3R chr4 chrX"

replicates="merge"

mcoolFileDir=./mcoolFiles/

n=0
for condition in ${conditions} ;
do
    for replicate in ${replicates} ;
    do
	for resolution in ${resolutions} ;
	do
	    for chrom1 in ${chroms} ;
	    do
		for maxDist in $maxDists ;
		do
		    chrLength=$(cat ./scripts/chrom_sizes.txt | awk -v c=${chrom1} '{if("chr"$1==c){print $2}}')	    
                    if [[ ${chrLength} -lt ${maxDist} ]];                                                          
                    then                                                                                           
			maxDist=${chrLength}                                                                           
			outDir=calledLoops_${resolution}_upTo_${maxDist}bp                                             
			mkdir -p ${outDir}
                    fi 
		    
		    outDir=calledLoops_${resolution}_upTo_${maxDist}bp
		    mkdir -p ${outDir}
		    for mcoolFile in $(ls -1 ${mcoolFileDir}/*${condition}_*${replicate}*.mcool) ;
		    do
			mcoolFileName=$(echo $mcoolFile | sed "s,/, ,g" | awk '{print $NF}')
			for sparsityThreshold in ${sparsityThresholds} ;
			do
			    for sigmaZero in $sigmaZeros ;		    
			    do
				outFile=${outDir}/${mcoolFileName%.mcool}_${chrom1}_${chrom1}_loops_${resolution}_sparsityThreshold_${sparsityThreshold}_sigmaZero_${sigmaZero}_iteration_5
				if [[ ! -e ${outFile} ]];
				then
				    n=$(bash ~/TOOLS/running_jobs.sh | tail -2 | awk '{s+=$NF}END{print s}')
				    while [[ ${n} -ge 2500 ]];
				    do					
					echo ${condition} ${resolution} ${sparsityThreshold} ${maxDist} ${sigmaZero} ${chrom1} ${n} $(date)
					sleep 60s
					n=$(bash ~/TOOLS/running_jobs.sh | tail -2 | awk '{s+=$NF}END{print s}')
				    done
				    
				    echo ${condition} ${resolution} ${sparsityThreshold} ${maxDist} ${sigmaZero} ${chrom1} ${n}
				    sbatch scripts/02_callLoops_using_mustache.sh ${condition} ${resolution} ${sparsityThreshold} ${maxDist} ${sigmaZero} ${chrom1}
				else
				    ls -lrtha ${outFile}
				fi
			    done
			done
		    done
		done
	    done
	done
    done
done
