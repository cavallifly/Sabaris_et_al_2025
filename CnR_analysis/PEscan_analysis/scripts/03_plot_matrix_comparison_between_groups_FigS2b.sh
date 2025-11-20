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

for microc in Embryo WD ED;
do
    echo $microc

    rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/*allPeaks_development_k_6_cl?*peScan_text.mtx .

    bash ../scripts/plot_matrix.sh    

    outFile=distPEscan_mC_${microc}_E_ED_WD_allPeaks_development_k_6_min${minDist}_max${maxDist}_range10_res250_peScan_centralRegion_8pixels_heatmaps_data.txt
    rm -fr ${outFile}
    for file in $(ls -1 *allPeaks_development_k_6_cl?_*_centralRegion_8pixels_heatmaps_data.txt);
    do
	echo $file
	cluster=$(echo $file | sed "s/_/ /g" | awk '{print $13}')
	echo $cluster
	cat ${file} | grep -w ${microc} | awk -v c=${cluster} '{print c,$2}' >> ${outFile}
    done
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n
done

Rscript ../scripts/03_violinPlots_within_clusters_scores_FigS2b.R 0

