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
    
for microc in Embryo WD ;
do
    echo $microc
    peak=WD_GAFs
    random=WD_random
    if [[ $microc == "Embryo" ]];
    then
	peak=E_GAFs
	random=EMBRYO_random
    fi

    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_mC_*FC2_${peak}_*peScan_text.mtx .
    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_mC_*q*_${peak}_*peScan_text.mtx .
    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/*_${random}_*peScan_text.mtx .

    bash ../scripts/plot_matrix.sh

    outFile=distPEscan_mC_${microc}_min${minDist}_max${maxDist}_range10_res250_peScan_centralRegion_8pixels_heatmaps_data.txt
    for file in $(ls -1 distPEscan_mC_*FC2_${peak}_*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cat ${file} | grep -w ${microc} | awk '{print "All",$2}' > ${outFile}
    done
    for file in $(ls -1 distPEscan_mC_*q*_${peak}_*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	quartile=$(echo $file | sed "s/_/ /g" | awk '{print $9}' | sed "s/q/Q/g")
	cat ${file} | grep -w ${microc} | awk -v q=${quartile} '{print q,$2}' >> ${outFile}
    done
    for file in $(ls -1 *_${random}_*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cat ${file} | grep -w ${microc} | awk '{print "Random",$2}' >> ${outFile}
    done
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n

done

conda run -n DEseq2 Rscript ../scripts/violinPlots_within_quartiles_scores_Fig1f.R 0
#Rscript ../scripts/violinPlots_within_quartiles_scores_Fig1f.R 0
