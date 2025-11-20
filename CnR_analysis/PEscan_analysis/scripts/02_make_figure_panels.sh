minDist=10
maxDist=500

mkdir -p Figures_minDist${minDist}_maxDist${maxDist}

for panel in Fig1cd FigS2b FigS3a FigS3b FigS3c ;
do
    
    outDir=material_for_${panel}_minDist${minDist}_maxDist${maxDist}
    if [[ -d ${outDir} ]];
    then
	echo "${outDir} exists. Remove it to repeat the analysis!"
	continue
    fi
    if [[ -d Figures_minDist${minDist}_maxDist${maxDist}/${outDir} ]];
    then
	echo "${outDir} exists in Figures_minDist${minDist}_maxDist${maxDist}. Remove it to repeat the analysis!"
	continue
    fi
    echo ${outDir}
    
    mkdir -p ${outDir}
    cd ${outDir}

    bash ../scripts/03_plot_matrix_comparison_between_groups_${panel}.sh ${minDist} ${maxDist} &> 03_plot_matrix_comparison_between_groups_${panel}.out
    
    ls -lrtha
    
    cd ..

    mv -v ${outDir} Figures_minDist${minDist}_maxDist${maxDist}
    echo

done

