rm -fvr _inFile_list_for_ANOVA

# Panel f
outFile=scoreMapsk250kexp500_r3000bp_NetA_Fig4and5.tsv
rm -fvr ${outFile}

files="scoreMapsk250kexp500_r3000bp_Fig5_larvae_DWT_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2Fab7_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26Pho_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2Virilis_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2x3end_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2en_pre2NetA_pre1NetA.tsv"

for file in $(ls -1 ${files}) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2NetA//g" -e "s/GAGA//g" -e "s/DPRE2_/DPRE2/g" -e "s/DPRE2Fab7_/Fab7/g" -e "s/DPRE2en_/en/g" -e "s/DPRE2hsp26Pho_/Hsp26+3xPHO/g" -e "s/DPRE2hsp26_/Hsp26/g" -e "s/DWT_/WT/g" -e "s/DPRE2Virilis_/Vir/g" -e "s/DPRE2x3end_/2xPRE1/g" -e "s/DWT_/WT/g" -e "s/DWT_/WT/g")    

    #wc -l $file
    echo $file $name
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

conda run -n DEseq2 Rscript scripts_clean/one-way_ANOVA_Fig4and5.R 0
