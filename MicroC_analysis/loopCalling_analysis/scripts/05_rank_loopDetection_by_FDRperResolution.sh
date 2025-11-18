
assembly=dm6
author=BeS

conditions="ED_PH18 ED_PH29 ED_PHD11 Embryo_WT WD_WT LD_WT"
resolutions="400 800 1000 2000 4000 8000 10000 20000 40000"
assay=microc

for condition in $conditions ;
do
    outFile=best_Nloops_per_FDR_threshold_${condition}_topNfraction.txt    
    outFileLoops=Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_PeriCentromere_FDR_resDependent_filtered_merged.bedpe
    rm -fr ${outFile} _tmp_${condition} ${outFileLoops}
    
    for resolution in $resolutions ;
    do
	sort -k 6,6n Nloops_per_FDR_threshold_${condition}.txt | awk -v r=$resolution '{if($2==r){print $0}}' | tail -1 >> ${outFile}
    done

    for resolution in $resolutions ;
    do
	FDR=$(awk -v r=$resolution '{if($2==r){print $3}}' ${outFile})

	awk '{print $NF}' loopFiles_per_FDR_${condition}/Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_res_${resolution}bp_PeriCentromere_FDR_${FDR}_filtered_merged_commonFurlong.bedpe >> _tmp_${condition}
	awk -v r=$resolution '{print $0,r}' loopFiles_per_FDR_${condition}/Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_res_${resolution}bp_PeriCentromere_FDR_${FDR}_filtered_merged.bedpe | grep -v "chrom1\|chrA" >> ${outFileLoops}
	
	sort _tmp_${condition} | uniq > _tmp_${condition}_1
	mv _tmp_${condition}_1 _tmp_${condition}
    done

    totLoops=$(awk '{s+=$4}END{print s}' ${outFile})
    sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n ${outFileLoops} > _tmp_${outFileLoops} ; mv _tmp_${outFileLoops} ${outFileLoops}
    wc -l ${outFileLoops}
    
    echo "Identified $totLoops loops and $(sort _tmp_${condition} | uniq | wc -l) Furlong loops over $(wc -l loopFurlongPlots_${condition}/common*/*bedpe)" >> ${outFile}
    cat $outFile
    echo
    mv _tmp_${condition} Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_PeriCentromere_FDR_resDependent_filtered_merged_commonFurlong.bedpe
    #cat Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_PeriCentromere_FDR_resDependent_filtered_merged_commonFurlong.bedpe
done
