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

rsync -avz ../distPEscan_mC_*_minDist${minDist}_maxDist${maxDist}/RData/distPEscan_mC_*AllE*renamed*ext_4000*peScan_text.mtx .

for file in $(ls -1 dist*AllE_AllE* 2> /dev/null);
do
    echo $file
    mv -v $file $(echo $file | sed "s/Loops_AllE_AllE_renamedLoops_//g")
done
ls -lrtha *peScan_text.mtx

bash ./scripts/03a_plot_matrix_overlappingLoopsBetweenTissues.sh

for file in $(ls -1 *centralRegion_${centralRegion}pixels*data.txt);
do
    category=$(echo $file | sed -e "s/_ext_/ /g" -e "s/000bp_/ /g" | awk '{print $3}')
    echo $file
    
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${file} | sed "s/PH18/ED/g" | sort -k 2,2n    
done

Rscript ./scripts/03b_quantification_overlappingLoopsBetweenTissues.R

cd ..
