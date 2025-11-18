ln -s ../../mustache_loopcaller/*default/Loop* .
outLog=merge_neighbour_loops_between_conditions.log
rm -fr _*

for resolution in 5000 10000 all ;
do

    # Merged loops    
    outFile=Loops_cis_chromosome_${resolution}bp_default_filtered_merged_merged.tsv
    if [[ ${resolution} == "all" ]];
    then
	outFile=Loops_cis_chromosome_default_filtered_merged_merged.tsv
    fi
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile}
	continue
    fi
    echo $outFile

    for chrom in chr2L chr2R chr3L chr3R chr4 chrX ;
    do
	echo $chrom >> ${outLog}
	
	# Get all the loops in a single file
	rm -fvr _${chrom}_loopFile   
	for loopFile in $(ls -1 Loops_scores_cis_chromosome_*${resolution}* 2> /dev/null);
	do
	    #echo $loopFile
	    cat ${loopFile} | grep ${chrom} | wc -l >> ${outLog}
	    cat ${loopFile} | grep ${chrom} >> _${chrom}_loopFile
	    #head ${loopFile}
	done # Close cycle over ${loopFile}	
	if [[ ${resolution} == "all" ]];
	then
	    for loopFile in $(ls -1 Loops_cis_chromosome_*bp* 2> /dev/null);
	    do
		#echo $loopFile
		cat ${loopFile} | grep ${chrom} | wc -l >> ${outLog}
		cat ${loopFile} | grep ${chrom} >> _${chrom}_loopFile
		#head ${loopFile}
	    done # Close cycle over ${loopFile}
	fi

	echo "Merging loops from different conditions overlapping or involving 1st neighbors bins!" >> ${outLog}
	# Sort the loops
	sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n _${chrom}_loopFile | uniq > _a ; mv _a _${chrom}_loopFile
	
	# Check if there is any loop to analyze
	nLoops=$(wc -l _${chrom}_loopFile | awk '{print $1}')
	echo "Initial loops $nLoops" >> ${outLog}	       
	if [[ $nLoops -eq 0 ]];
	then
	    rm -fvr _${chrom}_loopFile _${chrom}_merged _${chrom}_intersected _tmp >> ${outLog}
	    continue
	fi
	
	#https://github.com/billgreenwald/pgltools/issues/3
	#https://crazyhottommy.blogspot.com/2017/12/merge-enhancer-promoter-interaction.html
	if [[ ${resolution} != "all" ]];
	then
	    cat _${chrom}_loopFile | cut -f1-6 | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' | python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/merge.py -stdInA -d ${resolution} > _${chrom}_merged
	else
	    cat _${chrom}_loopFile | cut -f1-6 | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' | python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/merge.py -stdInA > _${chrom}_merged
	fi
	cat _${chrom}_merged | sort -k 1,1 -k 2,2n -k 4,4 -k 5,5n > _a ; mv _a _${chrom}_merged
	echo "Remaining loops ${chrom} $(wc -l _${chrom}_merged | awk '{print $1}')" >> ${outLog}
	
	cat _${chrom}_merged >> ${outFile}
	rm -fvr _${chrom}_merged _tmp >> ${outLog}
	echo >> ${outLog}
	echo "### DONE ###"
	
    done # Close cycle over $chrom       
done # Close cycle over ${resolution}


echo "Intersecting loops from different resolutions overlapping or involving 1st neighbors bins!" >> ${outLog}
# Intersected loops
outFileI=Loops_cis_chromosome_default_filtered_merged_intersected.tsv
if [[ -e ${outFileI} ]];
then
    ls -lrtha ${outFileI}
else
    echo $outFileI
    
    for chrom in chr2L chr2R chr3L chr3R chr4 chrX ;
    do
	echo $chrom >> ${outLog}
	
	# Get all the loops in two files A = 5kb and B = 10kb
	rm -fvr _${chrom}_loopFileA _${chrom}_loopFileB
	for loopFile in $(ls -1 Loops_scores_cis_chromosome_*5000bp* 2> /dev/null);
	do
	    #echo $loopFile
	    cat ${loopFile} | grep ${chrom} | wc -l >> ${outLog}
	    cat ${loopFile} | grep ${chrom} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n >> _${chrom}_loopFileA
	    #head ${loopFile}
	done # Close cycle over ${loopFile}
	for loopFile in $(ls -1 Loops_scores_cis_chromosome_*10000bp* 2> /dev/null);
	do
	    #echo $loopFile
	    cat ${loopFile} | grep ${chrom} | wc -l >> ${outLog}
	    cat ${loopFile} | grep ${chrom} | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n >> _${chrom}_loopFileB
	    #head ${loopFile}
	done # Close cycle over ${loopFile}
	
	#https://github.com/billgreenwald/pgltools/issues/3
	#https://crazyhottommy.blogspot.com/2017/12/merge-enhancer-promoter-interaction.html
	#python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/intersect.py -a _${chrom}_loopFileA -b _${chrom}_loopFileB > _${chrom}_intersected
	python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/intersect.py -a _${chrom}_loopFileB -b _${chrom}_loopFileA > _${chrom}_intersected	
	
	echo "Remaining loops ${chrom} $(wc -l _${chrom}_intersected | awk '{print $1}')" >> ${outLog}
	
	cat _${chrom}_intersected >> ${outFileI}
	rm -fvr _${chrom}_loopFileA _${chrom}_loopFileB _${chrom}_intersected _tmp >> ${outLog}
	echo >> ${outLog}	    
	
    done # Close cycle over $chrom       
fi
sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n ${outFileI} > _a ; mv _a ${outFileI}

echo "Intersecting loops from from our and other analysis!" >> ${outLog}
for refFile in ../../../2023_12_28_MolCell_Pollex_Furlong_A_high-confidence_set_of_prominent_chromatin_loops_mmc2.txt ;
do
    ls -lrtha ${refFile}
    #sort -k 2,2 -k 3,3n -k 4,4n ${refFile} | head
    extension=2500
    awk -v e=${extension} '{print $2,$3-e,$3+e,$2,$4-e,$4+e}' ${refFile} | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | grep -v ID | grep chr4 > _loopFileA
    head _loopFileA

    #for loopFile in ${outFileI} ${outFile} ;
    for loopFile in $(ls -1 ../../mustache_*/mC*/*tsv) #$(ls -1 Loop*);
    do
	echo $loopFile
	sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n ${loopFile} | grep -vi bin > _loopFileB
	head -30 _loopFileB

	wc -l _loopFileA _loopFileB
	python /home/michael.szalay/anaconda3/envs/mustache/lib/python3.8/site-packages/PyGLtools/intersect.py -a _loopFileA -b _loopFileB | wc -l
	
    done # Close cycle over ${loopFile}
done # Close cycle over ${refFile}    

