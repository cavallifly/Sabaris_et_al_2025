centralRegion=8


for filePH18 in $(ls -1 *mC_ED_PH18*range4*text.mtx);
do
    
    fileEmbryo=$( echo $filePH18 | sed "s/mC_ED_PH18/mC_Embryo/g")
    fileWD=$( echo $filePH18 | sed "s/mC_ED_PH18/mC_WD/g")    
    
    outFile=$(echo ${filePH18%_text.mtx}_centralRegion_${centralRegion}pixels_heatmaps.pdf | sed -e "s/ED_PH18_//g" -e "s/AllConditions_AllConditions_renamedLoops_//g")
    if [[ -e ${outFile} ]];
    then
	continue
    fi
    ls -lrtha ${filePH18} ${filePH29}
    
    echo

    awk '{for(i=2;i<=NF;i++){print NR-1,i-2,$i}}' <( awk '{if(NR>1) print $0}' ${filePH18})   > _tmp_matrix_PH18
    awk '{for(i=2;i<=NF;i++){print NR-1,i-2,$i}}' <( awk '{if(NR>1) print $0}' ${fileEmbryo}) > _tmp_matrix_Embryo
    awk '{for(i=2;i<=NF;i++){print NR-1,i-2,$i}}' <( awk '{if(NR>1) print $0}' ${fileWD})     > _tmp_matrix_WD
    
    head -1 _tmp_matrix_PH18
    tail -1 _tmp_matrix_PH18    
    size=$(tail -1 _tmp_matrix_PH18 | awk '{print $1}')
    echo $size
    
    cbmin=$(awk 'BEGIN{min=$3}{if($3<max){min=$3}}END{printf("%.2f",min)}' _tmp_matrix_*)        
    cbmax=$(awk 'BEGIN{max=$3}{if($3>max){max=$3}}END{printf("%.2f",max)}' _tmp_matrix_*)
    echo ${cbmin} ${cbmax}

    cbmin=-1.5 #$(echo ${cbmin} ${cbmax} | awk '{max=sqrt($1*$1); if(sqrt($2*$2)>max){max=sqrt($2*$2)}}END{print -max}')
    cbmax=1.5 #$(echo ${cbmin} ${cbmax} | awk '{max=sqrt($1*$1); if(sqrt($2*$2)>max){max=sqrt($2*$2)}}END{print +max}')    
    echo ${cbmin} ${cbmax}    

    name=$(echo ${filePH18} | sed -e "s/_PH18//g" -e "s/_res500_peScan_text\.mtx//g" -e "s/_res250_peScan_text\.mtx//g")
    echo ${filePH18}

    Nloops=$(echo ${filePH18} | sed -e "s/N1/ 1/g" -e "s/N2/ 2/g" -e "s/N3/ 3/g" -e "s/N4/ 4/g" -e "s/N5/ 5/g" -e "s/N6/ 6/g" -e "s/N7/ 7/g" -e "s/N8/ 8/g" -e "s/N9/ 9/g" | awk '{print $2}' | sed "s/_/ /g" | awk '{print $1}')
    echo "Number of loops "$Nloops
    
    xmin=$(echo ${filePH18} | sed -e "s/range/ /g" | awk '{print $2}' | sed "s/_/ /g" | awk '{print "-"$1"kb"}')
    xmax=$(echo ${filePH18} | sed -e "s/range/ /g" | awk '{print $2}' | sed "s/_/ /g" | awk '{print "+"$1"kb"}')
    ymin=$(echo ${filePH18} | sed -e "s/range/ /g" | awk '{print $2}' | sed "s/_/ /g" | awk '{print "-"$1"kb"}')
    ymax=$(echo ${filePH18} | sed -e "s/range/ /g" | awk '{print $2}' | sed "s/_/ /g" | awk '{print "+"$1"kb"}')        
    echo ${xmin} ${xmax} ${ymin} ${ymax}

    low=$( awk -v s=$size -v c=$centralRegion 'BEGIN{low=(s+1)/2-int(c/2); high=(s+1)/2+int(c/2); print low}')
    high=$(awk -v s=$size -v c=$centralRegion 'BEGIN{low=(s+1)/2-int(c/2); high=(s+1)/2+int(c/2); print high}')
    echo $low $high

    sed -e "s/XXXhighXXX/${high}/g" -e "s/XXXlowXXX/${low}/g" -e "s/XXXcbminXXX/${cbmin}/g"  -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXxminXXX/${xmin}/g" -e "s/XXXxmaxXXX/${xmax}/g" -e "s/XXXyminXXX/${ymin}/g" -e "s/XXXymaxXXX/${ymax}/g" -e "s/XXXsizeXXX/${size}/g" -e "s/XXXNXXX/${Nloops}/g" ./scripts/03a_plot_contact_matrix_overlappingLoopsBetweenTissues.gp | gnuplot

    bash ./scripts/ps2pdf.sh
    mv heatmap_map.pdf ${outFile}
    
    # Quantify differences
    rm -fvr _tmp.txt
    wc -l _tmp_matrix_*

    awk -v s=$size -v c=$centralRegion '{low=(s+1)/2-int(c/2); high=(s+1)/2+int(c/2); if((low<($1+1) && ($1+1)<=high) && (low<($2+1) && ($2+1)<=high)){print "PH18",$3}}'  _tmp_matrix_PH18 > _tmp.txt
    awk -v s=$size -v c=$centralRegion '{low=(s+1)/2-int(c/2); high=(s+1)/2+int(c/2); if((low<($1+1) && ($1+1)<=high) && (low<($2+1) && ($2+1)<=high)){print "Embryo",$3}}'  _tmp_matrix_Embryo >> _tmp.txt
    awk -v s=$size -v c=$centralRegion '{low=(s+1)/2-int(c/2); high=(s+1)/2+int(c/2); if((low<($1+1) && ($1+1)<=high) && (low<($2+1) && ($2+1)<=high)){print "WD",$3}}'  _tmp_matrix_WD >> _tmp.txt    
    wc -l _tmp.txt
    cp _tmp.txt ${outFile%.pdf}_data.txt 

    cbmin=$(awk 'BEGIN{min=$2}{if($2<max){min=$2}}END{printf("%.2f",min)}' _tmp.txt)
    echo $cbmin $(awk -v c=${cbmin} 'BEGIN{print c*1.05}')
done
