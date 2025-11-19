rm -fvr _inFile_list_for_ANOVA
for cellType in LD ED Embryo ;
do
    if [[ ${cellType} == "Embryo" ]];
    then
	cellType="Fig2_"${cellType}
    fi
    if [[ ${cellType} == "ED" ]];
    then
	cellType="Fig4_"${cellType}
    fi    
    if [[ ${cellType} == "LD" ]];
    then
	cellType="Fig5_"${cellType}
    fi    
    
    for PREloop in pre1dac pre1en pre1NetA ;
    do
	PREloopName=$(echo $PREloop | sed "s/1//g")
	echo $PREloop $PREloopName

	outFile=scoreMapsk250kexp500_r3000bp_${cellType}_${PREloopName}.tsv
	rm -fvr $outFile
	echo ${outFile} >> _inFile_list_for_ANOVA
	
	for file in $(ls -1 *${cellType}*${PREloop}*tsv);
	do
	    condition=$(echo $file | sed -e "s/_/ /g" | awk '{print $5}')
	    echo $file $condition
	    
	    awk -v c=${condition} '{print c,$NF}' ${file} | grep -v score >> ${outFile}
	    
	done
	echo
    done
done

for inFile in $(cat _inFile_list_for_ANOVA);
do

    sed "s/XXXinFileXXX/${inFile}/g" ./scripts_clean/one-way_ANOVA.R > _tmp.R

    #conda run -n DEseq2 Rscript _tmp.R
    Rscript _tmp.R
    rm -fvr _tmp.R

done
