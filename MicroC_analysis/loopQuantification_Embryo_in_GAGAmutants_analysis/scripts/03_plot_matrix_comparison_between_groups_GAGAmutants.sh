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

outFile=distPEscan_min${minDist}_max${maxDist}_range4_res250_peScan_centralRegion_${centralRegion}pixels_heatmaps_data.txt
rm -fvr $outFile
rsync -avz ../distPEscan_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_*peScan_text.mtx .
ls -lrtha distPEscan_*peScan_text.mtx

bash ./scripts/03a_plot_matrix_GAGAmutants.sh

file=$(ls -1 distPEscan_*_centralRegion_${centralRegion}pixels_heatmaps_data.txt)
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${file} | sort -k 1,1n


Rscript ./scripts/03b_quantification_overlappingLoopsBetweenTissues.R

