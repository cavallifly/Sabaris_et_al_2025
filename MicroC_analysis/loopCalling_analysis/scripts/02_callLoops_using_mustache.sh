#!/bin/bash

#SBATCH --job-name mustache
#SBATCH -n 10                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 150Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o mustache_500bp.out # File to which STDOUT will be written
#SBATCH -e mustache_500bp.out # File to which STDERR will be written 



# Locally, that is within 500000 bp separation between the loop anchors we can call loops
# using Hi-C maps of 1 and 2 kb resolution

# At larger separations, that is within 10000000 bp separation between the loop anchors we can call loops
# using Hi-C maps of 100, 20, 15, 10 and 5 kb resolution

#conditions="Embryo_WT WD_WT LD_WT ED_PH29 ED_PHD11 ED_PH18" # ED_PH29DStoPH29 ED_PHD11DStoPH29 ED_PH18DStoPH29"
conditions=$1
#resolutions="500 1000 2000 4000 8000"
resolutions=$2
#sparsityThresholds=$(seq 1.00 0.01 1.00) # <- sparsityThreshold == 1.00 gives very clean loops!
sparsityThresholds=$3
#maxDists="50000 100000 200000 500000 $(seq 1000000 1000000 10000000) $(seq 10000000 5000000 35000000)"
maxDist=$4
#sigmaZeros=$(seq 2.6 -0.1 0.6)
sigmaZeros=$5
#chroms="chr2L chr2R chr3L chr3R chrX"
chroms=$6

#sbatch scripts/01_callLoops_using_mustache.sh ${condition} ${resolution} ${sparsityThreshold} ${maxDist} ${sigmaZero} ${chrom1}

replicates="merge"

mcoolFileDir=/work/user/mdistefano/2022_08_22_microC_eyeDiscs_Bernd/01_cool_files/

for condition in ${conditions}
do
    for replicate in ${replicates}
    do
	for resolution in ${resolutions}
	do
	    #if [[ $resolution -eq 500 ]];
	    #then
	    #maxDist=$4
	    #fi

	    outDir=calledLoops_${resolution}_upTo_${maxDist}bp
	    mkdir -p ${outDir}

	    for mcoolFile in $(ls -1 ${mcoolFileDir}/*${condition}_*${replicate}*.mcool)
	    do
		mcoolFileName=$(echo $mcoolFile | sed "s,/, ,g" | awk '{print $NF}')
		echo $mcoolFile
		for sparsityThreshold in $sparsityThresholds ;
		do
		    #for sigmaZero in $(seq 0.6 0.1 2.6);
		    for sigmaZero in $sigmaZeros ;		    
		    do
			for iteration in 5 ; #$(seq 5 1 15); <- Iteration 5 10 15 gives the same results
			do	
			    for chrom1 in ${chroms} ;
			    do
				for chrom2 in ${chroms} ;
				do
				    chrLength=$(cat /work/user/mdistefano/mishaDB/trackdb/dm6/chrom_sizes.txt| awk -v c=${chrom1} '{if("chr"$1==c){print $2}}')
                                    #if [[ ${chrLength} -lt ${maxDist} ]];
                                    #then
				    #maxDist=${chrLength}
				    #outDir=calledLoops_${resolution}_upTo_${maxDist}bp
				    #mkdir -p ${outDir}
                                    #fi
				    
				    outFile=${outDir}/${mcoolFileName%.mcool}_${chrom1}_${chrom2}_loops_${resolution}_sparsityThreshold_${sparsityThreshold}_sigmaZero_${sigmaZero}_iteration_${iteration}
				    if [[ $chrom1 != $chrom2 ]];
				    then
					continue
				    fi
				    
				    if [[ -e ${outFile} ]];
				    then
					ls -lrtha ${outFile}
					continue
				    fi
				    touch ${outFile}
				    
				    echo "${resolution} upTo ${maxDist} bp"
				    python ./scripts/mustache/mustache/mustache.py -f ${mcoolFile} -o ${outFile} -r $resolution --processes 10 --octaves 2 -d ${maxDist} -ch ${chrom1} -ch2 ${chrom2} --pThreshold 0.1 --sparsityThreshold ${sparsityThreshold} --iteration ${iteration} --sigmaZero ${sigmaZero} &> ${outFile}.out
				    #cat mustacheRes_${outFile} >> ${outFile}.out
				    #mv mustacheRes_${outFile} ${outFile}
				    #exit
				done
			    done
			done
		    done
		done
	    done
	done
    done
done
