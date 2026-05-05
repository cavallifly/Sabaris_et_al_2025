minDist=3
maxDist=10000

mkdir -p Figures_minDist${minDist}_maxDist${maxDist}

for panel in overlappingLoopsBetweenEmbryos ;
do
    
    outDir=${panel}_minDist${minDist}_maxDist${maxDist}

    echo ${outDir}
    
    mkdir -p ${outDir}
    cd ${outDir}
    
    bash ../../scripts/03_plot_matrix_comparison_between_groups_${panel}.sh ${minDist} ${maxDist} &> 03_plot_matrix_comparison_between_groups_${panel}.out
    cd ..

    mv -v ${outDir} Figures_minDist${minDist}_maxDist${maxDist}
    echo

done

