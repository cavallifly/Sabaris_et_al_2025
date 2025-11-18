resolutions="400 800 1000 2000 4000 8000 10000 20000 40000"

for inFile in $(ls -1 Nloops_per_FDR_threshold_*.txt);
do
    condition=$(echo $inFile | sed -e "s/_/ /g" -e "s/\.txt//g" | awk '{print $(NF-1)"_"$NF}')
    echo $inFile $condition

    sort -k 2,2n -k 3,3n ${inFile} | uniq > _tmp ; mv _tmp ${inFile}
    
    for resolution in ${resolutions};
    do
	echo $resolution $(cat ${inFile} | awk -v r=$resolution '{if($2==r){h[$2]++}}END{for(i in h){print i,h[i]}}')
	cat ${inFile} | awk -v r=$resolution '{if($2==r){print $0}}' | awk '{h[$3]++}END{for(i in h){if(h[i]>1)print i,h[i]}}' | sort -k 1,1n > _tmp
	cat _tmp

	for FDR in $(awk '{print $1}' _tmp);
	do
	    echo $FDR $resolution
	    wc -l ${inFile}
	    awk -v FDR=$FDR -v r=$resolution '{if($2==r && $3==FDR){print $0}}' ${inFile}
	    awk -v FDR=$FDR -v r=$resolution '{if($2==r && $3==FDR){next}; print $0}' ${inFile} > _tmp_removed
	    mv _tmp_removed ${inFile}
	done

	for FDR in $(seq 0.0001 0.0001 0.1000);
	do
	    flag=$(awk -v FDR=$FDR -v r=$resolution 'BEGIN{c=0}{if($2==r && $3==FDR){c=1}}END{print c}' $inFile)

	    if [[ $flag -eq 0 ]];
	    then
		rm -fvr loopFiles_per_FDR_${condition}/Loops_*_contacts_*_merge_all_dm6_BeS_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_res_${resolution}bp_PeriCentromere_FDR_${FDR}_*.bedpe
	    fi
	done
	rm -fr _tmp
    done
done
