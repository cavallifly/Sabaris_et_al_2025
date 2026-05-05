# 0) params
#inFiles="Loops_Embryo_WT_BeS.bedpe Loops_Embryo_WT_Batut2022.bedpe Loops_Embryo_WT_Levo2022.bedpe Loops_Embryo_WT_Pollex2024.bedpe Loops_Embryo_WT_Dolsten2025.bedpe"
#set1=BeS
#set2=Dolsten2025

#set1=BeS
#set2=Pollex2024

#set1=BeS
#set2=Batut2022

#set1=Pollex2024
#set2=Dolsten2025

#set1=Batut2022
#set2=Dolsten2025

#rsync -avz ../Loops_Embryo_WT_${set1}.bedpe ../Loops_Embryo_WT_${set2}.bedpe .
#inFiles="Loops_Embryo_WT_${set1}.bedpe Loops_Embryo_WT_${set2}.bedpe"

set1=BeS
set2=Dolsten2025
set3=Pollex2024
set4=Batut2022

rsync -avz ../Loops_Embryo_WT_${set1}.bedpe ../Loops_Embryo_WT_${set2}.bedpe ../Loops_Embryo_WT_${set3}.bedpe ../Loops_Embryo_WT_${set4}.bedpe .
inFiles="Loops_Embryo_WT_${set1}.bedpe Loops_Embryo_WT_${set2}.bedpe Loops_Embryo_WT_${set3}.bedpe Loops_Embryo_WT_${set4}.bedpe"

echo $inFiles

for SLOP in 4000 ;
do
    rm -fvr _loops_*
    
    if [[ -e upset_plot_loop_overlap_ext_${SLOP}bp.pdf ]];
    then
	continue
    fi
    
    n=0
    rm -fr _all_loops _tmp_*
    for file in ${inFiles};
    do
	n=$((${n}+1))
	#cat $file | awk -v n=$n '{print $1,$2,$3,$4,$5,$6,"loopS"n"_"NR}' > _tmp_${file}
	cat $file | awk -v n=$n '{print $1,$2,$3,$4,$5,$6}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > _tmp_${file}
	cat _tmp_${file} | cut -f1-6 | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' | python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/merge.py -stdInA -d $((2*${SLOP})) | grep -v "#" | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | awk -v n=$n '{print $1,$2,$3,$4,$5,$6}' > _tmp_${file}_merged
	cat _tmp_${file} _tmp_${file}_merged | awk '{h[$0]++}END{for(i in h){if(h[i]==1){print i}}}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > ${file%.bedpe}_ext_${SLOP}bp_merged.bedpe

	awk '{print $0,sqrt(($2-$3)*($2-$3))}' ${file%.bedpe}_ext_${SLOP}bp_merged.bedpe | sort -k 7,7n | tac | awk '{if($7>1){n++; chrom1[n]=$1; start1[n]=$2; end1[n]=$3; chrom2[n]=$4; start2[n]=$5; end2[n]=$6;}else{for(i=0;i<=n;i++){if(chrom1[i]==$1 && chrom2[i]==$4){if((start1[i]<=$2 && $2<=end1[i]) && (start2[i]<=$5 && $5<=end2[i])){print $1,$2,$3,$4,$5,$6,"in",chrom1[i],int((start1[i]+end1[i])/2),int((start1[i]+end1[i])/2)+1,chrom2[i],int((start2[i]+end2[i])/2),int((start2[i]+end2[i])/2)+1}}}}}' > _tmp_${file%.bedpe}_ext_${SLOP}bp_merged.bedpe ; mv _tmp_${file%.bedpe}_ext_${SLOP}bp_merged.bedpe ${file%.bedpe}_ext_${SLOP}bp_merged.bedpe
	
	cat _tmp_${file}_merged | awk -v n=$n '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,$4,int(($5+$6)/2),int(($5+$6)/2)+1,"loopS"n"_"NR}' > _tmp_${file}
	rm -fvr _tmp_${file}_merged
    done

    cat _tmp_* >> _all_loops
    wc -l _tmp_* _all_loops
    #continue
    nA=0
    for A in $(ls -1 _tmp_*);
    do
	nA=$(($nA+1))
	nB=0
	for B in $(ls -1 _tmp_*);
	do
	    nB=$(($nB+1))
	    if [[ $nA -le $nB ]];
	    then
		continue
	    fi	    
	    OUTDIR=$(echo comparison_${A%.bedpe}_vs_${B%.bedpe}_overlap_ext_${SLOP} | sed "s/_tmp_//g")
	    mkdir -p $OUTDIR

	    # optional: how much to expand each anchor (in bp) for fuzzy matching
	    #SLOP=2000       # set to 0 to require direct overlap; adjust as needed
	    GENOME=$(ls -1 ../chrom_siz*)   # required only if SLOP>0; one-line per chrom: chr\tlength
	    
	    # 1) make anchor BEDs with loop ID and anchor-number
	    # anchor lines: chrom, start, end, loopID, anchor_num
	    #paste $A $B
	    awk 'BEGIN{OFS="\t"}{if($2<=$5){print $1,$2,$3,$7,1; print $4,$5,$6,$7,2;}else{print $4,$5,$6,$7,1; print $1,$2,$3,$7,2}}' $A > $OUTDIR/A_anchors.bed	    
	    #cat $OUTDIR/A_anchor1.bed $OUTDIR/A_anchor2.bed > $OUTDIR/A_anchors.bed
	    
	    awk 'BEGIN{OFS="\t"}{if($2<=$5){print $1,$2,$3,$7,1; print $4,$5,$6,$7,2;}else{print $4,$5,$6,$7,1; print $1,$2,$3,$7,2}}' $B > $OUTDIR/B_anchors.bed	    	    
	    #cat $OUTDIR/B_anchor1.bed $OUTDIR/B_anchor2.bed > $OUTDIR/B_anchors.bed
	    #paste $OUTDIR/B_anchors.bed $OUTDIR/A_anchors.bed
	    #exit
	    
	    # 2) optional: expand anchors by SLOP bp so slightly-offset anchors still match
	    if [[ ${SLOP} -gt 0 ]];
	    then
		head $OUTDIR/A_anchors.bed 
		bedtools slop -i $OUTDIR/A_anchors.bed -g $GENOME -b $SLOP > $OUTDIR/A_anchors.slop.bed
		bedtools slop -i $OUTDIR/B_anchors.bed -g $GENOME -b $SLOP > $OUTDIR/B_anchors.slop.bed
		A_anchors=$OUTDIR/A_anchors.slop.bed		
		B_anchors=$OUTDIR/B_anchors.slop.bed
		head $OUTDIR/A_anchors.slop.bed $OUTDIR/B_anchors.slop.bed
	    else
		A_anchors=$OUTDIR/A_anchors.bed
		B_anchors=$OUTDIR/B_anchors.bed
	    fi
	    echo $A_anchors $B_anchors
	    #paste $A_anchors $B_anchors
	    #exit
	    
	    # 3) intersect anchors (report both sides)
	    # Output columns: a_chr a_start a_end a_loopID a_anchorNum  b_chr b_start b_end b_loopID b_anchorNum
	    awk '{if($5==1){print $0}}' $A_anchors > _A_anchor1
	    awk '{if($5==1){print $0}}' $B_anchors > _B_anchor1
	    head _?_anchor1
	    bedtools intersect -a _A_anchor1 -b _B_anchor1 -wa -wb >  $OUTDIR/anchors_intersect.tsv
	    head $OUTDIR/anchors_intersect.tsv
	    
	    awk '{if($5==2){print $0}}' $A_anchors > _A_anchor2
	    awk '{if($5==2){print $0}}' $B_anchors > _B_anchor2
	    head _?_anchor2
	    bedtools intersect -a _A_anchor2 -b _B_anchor2 -wa -wb >> $OUTDIR/anchors_intersect.tsv	    
	    tail $OUTDIR/anchors_intersect.tsv
	    #exit

	    echo
	    # 4) aggregate: count anchors matched per (loopA, loopB) pair
	    # For each intersect line, extract loopA and loopB IDs and anchor numbers.
	    # Then count UNIQUE anchors per pair (in case of multi overlaps).
	    awk 'BEGIN{OFS="\t"}{if($5==$10){print $4,$9,$5,$10}}' $OUTDIR/anchors_intersect.tsv | sort -k1,1 -k2,2 | awk '{count[$1"_"$2]++}END{for(k in count) print k, count[k]}' | sort -k 2,2n > $OUTDIR/loopA_loopB_anchorCount.tsv
	    head $OUTDIR/loopA_loopB_anchorCount.tsv

	    echo
	    # 5) select pairs where count >= 2  (both anchors matched)
	    awk '{if($2==2){print $0}}' $OUTDIR/loopA_loopB_anchorCount.tsv | sort | uniq > $OUTDIR/matched_loop_pairs.tsv
	    head $OUTDIR/matched_loop_pairs.tsv
	    #exit
	    
	    # 6) summary stats	    
	    A_total=$(wc -l $A | awk '{print $1}')   # loops in A
	    B_total=$(wc -l $B | awk '{print $1}')   # loops in B	    
	    matched_pairs=$(wc -l $OUTDIR/matched_loop_pairs.tsv | awk '{print $1}')
	    A_matched=$(sed "s/_/ /g" $OUTDIR/matched_loop_pairs.tsv | sort -u | wc -l)
	    B_matched=$(sed "s/_/ /g" $OUTDIR/matched_loop_pairs.tsv | sort -u | wc -l)
	    
	    echo "loopsA_total: $A_total"
	    echo "loopsB_total: $B_total"
	    echo "matched_pairs (A <- B full-anchor matches): $matched_pairs"
	    echo "unique loops in A matched to something in B: $A_matched (${A_matched} / ${A_total} = $(echo ${A_matched} ${A_total} | awk '{printf("%.2f",$1/$2*100.)}'))"
	    echo "unique loops in B matched to something in A: $B_matched (${B_matched} / ${B_total} = $(echo ${B_matched} ${B_total} | awk '{printf("%.2f",$1/$2*100.)}'))"
	    #cat $OUTDIR/matched_loop_pairs.tsv
	    #exit
	    
	    # Update loop lists
	    cond1=$(head -1 $OUTDIR/matched_loop_pairs.tsv | sed -e "s/_/ /g" -e "s/loop//g" | awk '{print $1; print $3}' | sort | head -1 | awk '{print $1}')
	    cond2=$(head -1 $OUTDIR/matched_loop_pairs.tsv | sed -e "s/_/ /g" -e "s/loop//g" | awk '{print $1; print $3}' | sort | tail -1 | awk '{print $1}')
	    echo $cond1 $cond2
	    nToRename=$(wc -l $OUTDIR/matched_loop_pairs.tsv | awk '{print $1}')
	    cp $OUTDIR/matched_loop_pairs.tsv _matched_loop_pairs.tsv
	    
	    rm -fvr _changed_loops
	    for n in $(seq 1 1 ${nToRename});
	    do
		#awk -v n=$n '{if(NR==n){print $0}}' _matched_loop_pairs.tsv
		loopA=$(awk -v n=$n '{if(NR==n){print $0}}' _matched_loop_pairs.tsv | sed -e "s/_/ /g" | awk '{print $1"_"$2}')
		loopB=$(awk -v n=$n '{if(NR==n){print $0}}' _matched_loop_pairs.tsv | sed -e "s/_/ /g" | awk '{print $3"_"$4}')
		#nameNew=$(echo $loopA $loopB | sed -e "s/loop//g" -e "s/_/ /g" | awk '{print $1; print $3}' | sed "s/S/ /g" | awk '{for(i=1;i<=NF;i++){print "S"$i}}' | sort | uniq | awk '{printf("%s",$1)}')
		echo $loopA $loopB | sed -e "s/loop//g" -e "s/_/ /g" | awk '{print $1; print $3}'
		nameNew=$(echo $loopA $loopB | sed -e "s/loop//g" -e "s/_/ /g" | awk '{print $1; print $3}' | sed "s/S/ /g" | awk '{for(i=1;i<=NF;i++){print "S"$i}}' | sort -k 1,1n | uniq | awk '{printf("%s",$1)}')
		echo ${nameNew}

		touch _loops_${nameNew}
		loopNEW=$(awk -v loopA=$loopA -v loopB=$loopB 'BEGIN{v=0}{if($1==loopA || $1==loopB){print v=$2}}END{print v}' _changed_loops | uniq | tail -1)

		if [[ ${loopNEW} -eq 0 ]];
		then
		    awk -v name=${nameNew} 'BEGIN{n=0}{print $0;n+=1}END{n+=1;print "loop"name"_"n}' _loops_${nameNew} > _tmp_loops_${nameNew} ; mv _tmp_loops_${nameNew} _loops_${nameNew}
		    loopNEW=$(tail -1 _loops_${nameNew})
		else
		    echo $loopNEW
		fi   
		echo "New Loop $loopA $loopB ${loopNEW}"

		#cat ${A} | sed "s/\<${loopA}\>/${loopNEW}/g" > _tmp ; diff _tmp ${A} ; mv _tmp ${A}
		#cat ${B} | sed "s/\<${loopB}\>/${loopNEW}/g" > _tmp ; diff _tmp ${B} ; mv _tmp ${B}
		for file in $(ls -1 _tmp_*)
		do
		    #if [[ ${file} == ${A} ]];
		    #then
		    #continue
		    #fi
		    #if [[ ${file} == ${B} ]];
		    #then
		    #		continue
		    #fi
		    cat ${file} | sed "s/\<${loopA}\>/${loopNEW}/g" > _tmp ; diff _tmp ${file} ; mv _tmp ${file}
		    cat ${file} | sed "s/\<${loopB}\>/${loopNEW}/g" > _tmp ; diff _tmp ${file} ; mv _tmp ${file}
		done
		cat _all_loops | sed -e "s/\<${loopA}\>/${loopNEW}/g" -e "s/\<${loopB}\>/${loopNEW}/g" > _tmp ; diff _tmp _all_loops ; mv _tmp _all_loops

		echo ${loopA} ${loopNEW} >> _changed_loops
		echo ${loopB} ${loopNEW} >> _changed_loops
		tail _changed_loops
		#exit
		
		#cp _matched_loop_pairs.tsv _matched_loop_pairs_old.tsv
		#sed "s/_l/ l/g" _matched_loop_pairs.tsv | sed -e "s/\<${loopA}\>/${loopNEW}/g" -e "s/\<${loopB}\>/${loopNEW}/g" | awk '{print $1"_"$2,$3}' > _tmp ; diff _tmp _matched_loop_pairs.tsv ; mv _tmp _matched_loop_pairs.tsv
		#echo "Change in _matched_loop_pairs_old.tsv"
		#diff _matched_loop_pairs.tsv _matched_loop_pairs_old.tsv		
		#exit
	    done
	    
	    # matched_loop_pairs.tsv has columns: loopA  loopB  anchorCount (>=2)
	    # You can join these back to original BEDPEs to export matched loop coordinates.
	done
    done
    #rm -fr __tmp*
    
    sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n _all_loops > $OUTDIR/all_renamedLoops_ext_${SLOP}bp.bedpe
    rm -fr _all_loops
    for file in $(ls -1 _tmp_*);
    do
	echo $file
	
	mergeFile=$(echo ${file} | sed -e "s/_tmp_//g" -e "s/\.bedpe/_ext_${SLOP}bp_merged\.bedpe/g")
	echo $mergeFile
	
	outFile=$OUTDIR/$(echo $file | sed -e "s/_tmp_//g" -e "s/\.bedpe/_renamedLoops_ext_${SLOP}bp\.bedpe/g")
	echo $outFile
	echo
	cat ${mergeFile} $file | awk '{if(NF>7){h[$8"_"$9"_"$10"_"$11"_"$12"_"$13]=1; c[$1"_"$2"_"$3"_"$4"_"$5"_"$6]=$8"_"$9"_"$10"_"$11"_"$12"_"$13}else{loop=$1"_"$2"_"$3"_"$4"_"$5"_"$6; if(loop in h){for(j in c){if(c[j]==loop){print j,$7}}}else{print $0}}}' | sed "s/_/ /g" | awk '{print $1,$2,$3,$4,$5,$6,$7"_"$8}' > _tmp
	diff _tmp $file

	cat ${file} | sort -k 1,1 -k 2,2n  -k 3,3n -k 4,4 -k 5,5n  -k 6,6n > ${outFile}
	wc -l ${file} ${outFile}
	cp ${outFile} $(echo $file | sed -e "s/_tmp_//g" -e "s/\.bedpe/_renamedLoops_ext_${SLOP}bp\.bedpe/g")
	rm -fr $file
    done

    Rscript ./scripts/01_make_upset_plot_embryo.R $SLOP
    mv upset_plot_loop_overlap.pdf upset_plot_loop_overlap_ext_${SLOP}bp.pdf

    #cat Loops_Embryo_WT_${set1}_renamedLoops_ext_4000bp.bedpe | grep    S1S2 > Loops_Embryo_WT_${set1}-${set2}_renamedLoops_ext_4000bp_in_${set1}.bedpe
    #cat Loops_Embryo_WT_${set1}_renamedLoops_ext_4000bp.bedpe | grep -v S1S2 > Loops_Embryo_WT_${set1}_renamedLoops_ext_4000bp_notIn_${set2}.bedpe
    #wc -l Loops_Embryo_WT_${set1}-${set2}_renamedLoops_ext_4000bp_in_${set1}.bedpe Loops_Embryo_WT_${set1}_renamedLoops_ext_4000bp_notIn_${set2}.bedpe
    
    #cat Loops_Embryo_WT_${set2}_renamedLoops_ext_4000bp.bedpe | grep    S1S2 > Loops_Embryo_WT_${set1}-${set2}_renamedLoops_ext_4000bp_in_${set2}.bedpe
    #cat Loops_Embryo_WT_${set2}_renamedLoops_ext_4000bp.bedpe | grep -v S1S2 > Loops_Embryo_WT_${set2}_renamedLoops_ext_4000bp_notIn_${set1}.bedpe        
    #wc -l Loops_Embryo_WT_${set1}-${set2}_renamedLoops_ext_4000bp_in_${set2}.bedpe Loops_Embryo_WT_${set2}_renamedLoops_ext_4000bp_notIn_${set1}.bedpe

    rm -fvr Rplots.pdf
    
done

rm -fvr _*
