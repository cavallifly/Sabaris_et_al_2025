for resolution in 4000 ;
do

    # Overlap among all conditions
    inFiles=$(ls -1 Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_BeS-Batut2020-Levo2022-Pollex2024-Dolsten2025.bedpe)
    outFile=Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_allSets_in_BeS.bedpe
    wc -l ${inFiles}
    echo "${inFiles} -> ${outFile}"
    awk '{if($NF=="BeS"){print $0}}' ${inFiles} > ${outFile}
    wc -l ${outFile}
    echo
    
    for set1 in BeS Batut2020 Levo2022 Pollex2024 ;
    do
	for set2 in BeS Dolsten2025 Pollex2024 Batut2020 Levo2022 ;
	do
	    check=$(ls -1 Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_${set2}-${set1}_in_${set2}.bedpe 2> /dev/null | wc -l)
	    if [[ ${check} -ne 0 ]];
	    then
		continue
	    fi
	    
	    if [[ ${set1} == ${set2} ]];
	    then
		continue
	    fi
	    # Overlap ${set1} and ${set2}
	    inFiles=$(ls -1 Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_*${set1}*${set2}*.bedpe Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_*${set2}*${set1}*.bedpe 2> /dev/null | grep -v _in_ | grep -v _notIn_)
	    outFile=Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_${set1}-${set2}_in_${set1}.bedpe
	    cat ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | grep "${set1}\|${set2}"
	    wc -l Loops_Embryo_WT_${set1}_renamedLoops_ext_${resolution}bp.bedpe Loops_Embryo_WT_${set2}_renamedLoops_ext_${resolution}bp.bedpe
	    echo "${inFiles} -> ${outFile}"
	    awk -v s=$set1 '{if($NF==s){print $0}}' ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq > ${outFile}
	    wc -l ${outFile}
	    echo

	    # ${set1} and not in ${set2}
	    inFiles=$(ls -1 Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_*${set1}*.bedpe | grep -v ${set2} | grep -v _in_ |  grep -v _notIn_)
	    outFile=Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_${set1}_notIn_${set2}.bedpe
	    cat ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | grep "${set1}"
	    wc -l Loops_Embryo_WT_${set1}_renamedLoops_ext_${resolution}bp.bedpe
	    echo "${inFiles} -> ${outFile}"
	    awk -v s=${set1} '{if($NF==s){print $0}}' ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq > ${outFile}
	    wc -l ${outFile}
	    echo
	    
	    # ${set2} and not in ${set1}
	    inFiles=$(ls -1 Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_*${set2}*.bedpe | grep -v ${set1} | grep -v _in_ | grep -v _notIn_)
	    outFile=Loops_AllE_AllE_renamedLoops_ext_${resolution}bp_${set2}_notIn_${set1}.bedpe
	    cat ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | grep "${set2}"
	    wc -l Loops_Embryo_WT_${set2}_renamedLoops_ext_${resolution}bp.bedpe
	    echo "${inFiles} -> ${outFile}"
	    awk -v s=${set2} '{if($NF==s){print $0}}' ${inFiles} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq > ${outFile}
	    wc -l ${outFile}
	    echo
	done # Close cycle over $set1
	echo
    done # Close cycle over $set1

done
