centralRegion=8

minDist=3 #$1
maxDist=10000 #$2
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

mkdir Figures_minDist${minDist}_maxDist${maxDist}
cd Figures_minDist${minDist}_maxDist${maxDist}

outFile=distPEscan_min${minDist}_max${maxDist}_range4_res250_peScan_centralRegion_8pixels_heatmaps_data.txt
rm -fvr $outFile
for microc in ED_PH18 ;
do
    echo $microc

    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_mC_**peScan_text.mtx .
    rm -fvr *renamed*
    
    conda run -n misha bash ../scripts/plot_matrix.sh

    for file in $(ls -1 distPEscan_mC_*ALLl*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cat ${file} | awk '{print "All_"$1,$2}'   | sed -e "s/ED_PH18/Control/g" -e "s/ED_PH29/PH_KD/g" >> ${outFile}
    done
    for file in $(ls -1 distPEscan_mC_*noP*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cat ${file} | awk '{print "noPcG_"$1,$2}' | sed -e "s/ED_PH18/Control/g" -e "s/ED_PH29/PH_KD/g" >> ${outFile}
    done
    for file in $(ls -1 distPEscan_mC_*_P*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cat ${file} | awk '{print "PcG_"$1,$2}'   | sed -e "s/ED_PH18/Control/g" -e "s/ED_PH29/PH_KD/g" >> ${outFile}
    done
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n

done

#conda run -n DEseq2 Rscript ../scripts/violinPlots_within_quartiles_scores_Fig3d.R 0
#Rscript ../scripts/violinPlots_within_clusters_scores_Fig3d.R 0
conda run -n DEseq2 Rscript ../scripts/one-way_ANOVA_Fig3d.R

cd ..
