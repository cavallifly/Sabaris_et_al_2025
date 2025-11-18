#!/bin/bash

#SBATCH --job-name 04_filtering
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 04_filtering_loops_by_FDRperResolution.out # File to which STDOUT will be written
#SBATCH -e 04_filtering_loops_by_FDRperResolution.out # File to which STDERR will be written

assembly=dm6
author=BeS

conditions="larvae_DWT"
resolutions="1000 2000 4000 8000 10000 20000 40000"
assay=hic

#conditions="ED_PH18 ED_PH29 ED_PHD11 Embryo_WT WD_WT LD_WT"
#conditions="ED_PH18 ED_PH29 Embryo_WT WD_WT"
#resolutions="400 800 1000 2000 4000 8000 10000 20000 40000"
#assay=microc


#for condition in $1 ;
for condition in ${conditions} ; 
do
    FurlongLoops=loopFurlongPlots_${condition}/commonLoops/Furlong_commonLoops_${condition}.bedpe
    nLoops=$(cat ${FurlongLoops} | wc -l | awk '{print $1}')
    echo "Found ${nLoops} target loops in ${condition}"
    mkdir -p loopFiles_per_FDR_${condition}
    
    #for resolution in $2 ;
    for resolution in ${resolutions} ;    		      
    do

	for FDR in $(seq 0.0001 0.0001 0.1) ;
	do
	    echo $FDR
	    outFile=./loopFiles_per_FDR_${condition}/Loops_${assay}_contacts_${condition}_merge_all_${assembly}_${author}_allChromosomes_loops_allResolutions_sparsityThreshold_1.00_sigmaZero_merge_iteration_5_res_${resolution}bp_PeriCentromere_FDR_${FDR}_filtered_merged.bedpe
	    if [[ ! -e ${outFile} ]];
	    then
		touch ${outFile}
		awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%f\n",$1,$2,$3,$4,$5,$6,$7)}' calledLoops_${resolution}_upTo_*/*${condition}*5_PeriCentromere_filtered | grep -v BIN1 | awk -v t=$FDR '{if($7<t){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq | python ~/PyGLtools/merge.py -d ${resolution} -stdInA 2> /dev/null | grep -v "chrom1\|chrA" | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq >> ${outFile}
	    else
		continue
	    fi
	    cat ${outFile} | grep -v "chrom1\|chrA" > ${outFile}_tmp
	    mv ${outFile}_tmp ${outFile}
	    cat ${outFile}
	    nLoops=$(grep -v "chrom1\|chrA" ${outFile} | wc -l | awk '{print $1}')
	    
	    if [[ ! -e ${outFile%.bedpe}_commonFurlong.bedpe ]];
	    then
		touch ${outFile%.bedpe}_commonFurlong.bedpe
		for loop in $(seq 1 1 $nLoops);
		do		    
		    chrom1=$(awk -v l=$loop '{if(NR==l){print $1}}' ${FurlongLoops})
		    start1=$(awk -v l=$loop '{if(NR==l){print $2}}' ${FurlongLoops})
		    end1=$(awk -v l=$loop   '{if(NR==l){print $3}}' ${FurlongLoops})
		    chrom2=$(awk -v l=$loop '{if(NR==l){print $4}}' ${FurlongLoops})
		    start2=$(awk -v l=$loop '{if(NR==l){print $5}}' ${FurlongLoops})
		    end2=$(awk -v l=$loop   '{if(NR==l){print $6}}' ${FurlongLoops})
		    name=$(awk -v l=$loop   '{if(NR==l){print $7}}' ${FurlongLoops})

		    awk -v l=$loop   '{if(NR==l){print $0}}' ${FurlongLoops}
		    
		    awk -v r=${resolution} -v c1=${chrom1} -v s1=${start1} -v e1=${end1} -v c2=${chrom2} -v s2=${start2} -v e2=${end2} -v n=${name} '{if(c1!=$1){next}; sumM=$2+$3+$5+$6; sumF=s1+e1+s2+e2; if(sqrt((sumF-sumM)*(sumF-sumM))>1000000){next}; for(i=s1;i<=e1;i+=100){for(j=s2;j<=e2;j+=100){if((($2-0*r)<=i && i<=($3+0*r)) && (($5-0*r)<=j && j<=($6+0*r))){print $0,n;exit}}}}' <(grep -w ${chrom1} ${outFile}) | grep -v "chrom1\|chrA" >> ${outFile%.bedpe}_commonFurlong.bedpe
		done
	    fi
	    sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n ${outFile%.bedpe}_commonFurlong.bedpe | grep -v "chrom1\|chrA" | uniq > ${outFile%.bedpe}_commonFurlong.bedpe_tmp
            mv ${outFile%.bedpe}_commonFurlong.bedpe_tmp ${outFile%.bedpe}_commonFurlong.bedpe
	    
	    echo $condition $resolution $FDR $nLoops $(wc -l ${outFile%.bedpe}_commonFurlong.bedpe | awk '{print $1}') | awk '{if($4==0){print $0,0}else{print $0,$5/$4}}' >> Nloops_per_FDR_threshold_${condition}.txt
	    #exit
	done
    done
done
