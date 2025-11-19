conditions=$1
resolutions=$2
sparsityThresholds=$3
maxDist=$4
sigmaZeros=$5
chroms=$6

replicates="merge"

mcoolFileDir=./mcoolFiles/

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
