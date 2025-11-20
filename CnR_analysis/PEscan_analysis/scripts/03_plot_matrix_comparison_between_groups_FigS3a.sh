centralRegion=8

minDist=$1
maxDist=$2
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

for microc in Embryo;
do
    echo $microc
    if [[ $microc == "Embryo" ]];
    then
	peak=E_GAFs
	random=EMBRYO_random
    fi

    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/*E_GAF_allPeaks_Histones_k_3_k?_*peScan_text.mtx .

    bash ../scripts/04_plot_heatmaps.sh
    
    outFile=distPEscan_mC_${microc}_E_GAF_allPeaks_Histones_k_3_min${minDist}_max${maxDist}_range10_res250_peScan_centralRegion_8pixels_heatmaps_data.txt
    rm -fr ${outFile}
    for file in $(ls -1 *_GAF_allPeaks_Histones_k_3_k?_*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cluster=$(echo $file | sed "s/_/ /g" | awk '{print $11}' | sed "s/k/cl/g")
	cat ${file} | grep -w ${microc} | awk -v c=${cluster} '{print c,$2}' >> ${outFile}
    done
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n

done
    
Rscript ../scripts/03_violinPlots_within_clusters_scores_FigS3.R 0
