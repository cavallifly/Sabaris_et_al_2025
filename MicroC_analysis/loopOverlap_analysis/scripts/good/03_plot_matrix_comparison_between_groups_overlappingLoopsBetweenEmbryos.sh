centralRegion=8

minDist=3
maxDist=10000
if [[ ${minDist} == "" ]];
then
    echo "Please, provide the minDist as the first in-line parameter"
    exit
fi
if [[ ${maxDist} == "" ]];
then
    echo "Please, provide the maxDist as the first in-line parameter"
    exit
fi


rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_mC_**renamed*ext_4000*peScan_text.mtx .

for file in $(ls -1 dist*ext_4000*peScan_text.mtx);
do
    echo $file
    mv -v $file $(echo $file | sed "s/_Loops_Embryo_WT//g")
done
ls -lrtha *peScan_text.mtx

bash ./scripts/03a_plot_matrix_overlappingLoopsBetweenEmbryos.sh

for set in Dolsten2025 Batut2022 Pollex2024 ;
do

    outFile=distPEscan_mC_Embryo_minDist3_maxDist10000_ext_4000bp_BeS-${set}_min3_max10000_range4_res250_centralRegion_${centralRegion}pixels_data.txt
    rm -fvr $outFile
    for file in $(ls -1 *BeS*${set}*centralRegion_${centralRegion}pixels*data.txt);
    do
	category=$(echo $file | sed "s/_/ /g" | awk '{print $6"_"$10"_"$11}')
	echo $file
	echo $category
	
	cat $file | grep Embryo | awk -v c=$category '{print c,$2}' | sed "s/BeS-${set}_in_BeS/overlappingLoops/g" >> ${outFile}    
    done
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 2,2n

done
Rscript ./scripts/03b_quantification_overlappingLoopsBetweenEmbryos.R

cd ..
#done
