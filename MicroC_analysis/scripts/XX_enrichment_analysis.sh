minDist=10000
for maxDistance in 1000000 2000000 5000000
do

    for chrom in chr2L chr2R chr3L chrX chr3R;
    do
	
	inFile=scores_trans1Dinterval_${chrom}_all_genes_vs_all_domains.tab
	if [[ ! -e _${inFile} ]];
	then
	    (
		awk '{if(NF==1){r[$1]=$1}else{if($7 in r || $8 in r){next}else{print $0}}}' domains_to_remove.txt ${inFile} > _${inFile}
		#diff ${inFile} _${inFile}
		wc -l $inFile _${inFile}
	    ) &
	fi
	head $inFile
	
	for state1 in PcG Null Active ;
	do
	    l1=$(echo $state1 | wc | awk '{print $NF-1}')
	    echo $state1 $l1
	    
	    for state2 in PcG Null Active ;
	    do
		l2=$(echo $state2 | wc | awk '{print $NF-1}')
		echo $state2 $l2
		
		minDist=10000
		for maxDist in $(seq 100000 100000 ${maxDistance});
		do
		    outFile=scores_trans1Dinterval_${chrom}_${state1}_vs_${state2}_minDist_${minDist}bp_maxDist_${maxDist}bp.tsv
		    
		    if [[ -e ${outFile} ]];
		    then
			minDist=$maxDist		
			continue
		    fi
		    
		    (
			
			awk -v minDist=$minDist -v maxDist=$maxDist -v l1=$l1 -v state1=$state1 -v l2=$l2 -v state2=$state2 '{if((substr($8,1,l1)==state1 && substr($7,1,l2)==state2) || (substr($7,1,l1)==state1 && substr($8,1,l2)==state2)){c1=int(($3+$2)/2); c2=int(($6+$5)/2); dc=c1-c2; d=sqrt(dc*dc); if(minDist<=d && d<=maxDist){print $0}}}' _${inFile} > ${outFile}
			
			wc -l ${outFile}
		    ) &
		    minDist=$maxDist
		    
		done # Close cycle over $maxDist
		if [[ ${state1} == ${state2} ]];
		then
		    break
		fi	       
	    done # Close cycle over $state2
	done # Close cycle over $state1
	wait
    done # Close cycle over $chrom
    #exit
    
    #rm -fvr scores_trans1Dinterval_*_data.tsv
    width=100000
    for comparison in PcG_vs_PcG Active_vs_Active Null_vs_Null ;    
    do
	echo $comparison
	for chrom in chr3R All; # chr2L chr2R chr3L chr3R chrX All;
	do
	    chromName=${chrom}
	    outFile=scores_trans1Dinterval_${chromName}_${comparison}_every_${width}_upTo_${maxDistance}_data.tsv	    
	    if [[ ${chrom} == "All" ]];
	    then
		chrom=""
	    fi
	    #if [[ -e $outFile ]];
	    #then
	    #    continue
	    #fi
	    
	    for file in $(ls -1 *${chrom}*${comparison}*minDi*tsv);
	    do
		echo ${file} ${chrom} ${chromName}
		for condition in PH18 PH29 ;
		do		
		    #range=$(echo $file | sed "s/_/ /g" | awk -v w=${width} '{ind1=int($8/w)*w; ind2=int($10/w)*w; if($8<w){ind1=0; ind2=w}; if(ind1==ind2){ind2=ind1+w}; print ind1"-"ind2}' | sed -e "s/bp//g" -e "s/\.tsv//g")
		    range=$(echo $file | sed "s/_/ /g" | awk -v w=${width} '{ind1=int($8/w)*w; ind2=int($10/w)*w; if(ind1==ind2){ind2=ind1+w}; print ind1"-"ind2}' | sed -e "s/bp//g" -e "s/\.tsv//g")				
		    echo $condition $range
		    cat $file | grep $condition | awk '{h[$7"-"$8]+=$9; cnt[$7"-"$8]++}END{for(i in h) print i,h[i]/cnt[i],cnt[i]}' | awk -v c=${condition} -v r=${range} '{print c"_"r,$2}' >> ${outFile}
		    
		    #cat $file | grep $condition | awk '{h[$7"-"$8]+=$9; cnt[$7"-"$8]++}END{for(i in h) print h[i]/cnt[i]}' | ~/howmuch.csh
		done
	    done
	done
    done
    #exit
    
    for thresholdDist in 500000 400000 300000 200000 100000 ;
    do
	for comparison in PcG_vs_PcG Active_vs_Active Null_vs_Null ;    
	do
	    echo $comparison
	    for chrom in chr3R All ; #chr2L chr2R chr3L chr3R chrX All;
	    do
		chromName=${chrom}
		outFile=scores_trans1Dinterval_${chromName}_${comparison}_thresholdDist_${thresholdDist}_upTo_${maxDistance}_data.tsv
		if [[ ${chrom} == "All" ]];
		then
		    chrom=""
		fi
		if [[ -e $outFile ]];
		then
		    continue
		fi
		
		for file in $(ls -1 *${chrom}*${comparison}*minDi*tsv);
		do
		    echo ${file} ${chrom} ${chromName}
		    for condition in PH18 PH29 ;
		    do
			range=$(echo $file | sed -e "s/_/ /g" -e "s/bp//g" -e "s/\.tsv//g" | awk -v m=${maxDistance} -v t=$thresholdDist '{if($10<=t){print "10000-"t}else{print t"-"m}}')	    
			echo $condition $range
			
			cat $file | grep $condition | awk '{h[$7"-"$8]+=$9; cnt[$7"-"$8]++}END{for(i in h) print i,h[i]/cnt[i],cnt[i]}' | awk -v c=${condition} -v r=${range} '{print c"_"r,$2}' >> ${outFile}
			
			#cat $file | grep $condition | awk '{h[$7"-"$8]+=$9; cnt[$7"-"$8]++}END{for(i in h) print h[i]/cnt[i]}' | ~/howmuch.csh
		    done
		done
	    done
	done
    done
done
