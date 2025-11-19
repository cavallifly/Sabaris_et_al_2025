rm -fvr _inFile_list_for_ANOVA

# Panel f
outFile=scoreMapsk250kexp500_r3000bp_dac_Fig4f.tsv
rm -fvr ${outFile}

files="scoreMapsk250kexp500_r3000bp_Fig5_larvae_DWT_pre2dac_pre1dac.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26Pho_pre2dac_pre1dac.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26_pre2dac_pre1dac.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2_pre2dac_pre1dac.tsv"

for file in $(ls -1 ${files}) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2dac//g" -e "s/GAGA//g" -e "s/DPRE2_/DPRE2/g" -e "s/DPRE2hsp26Pho_/Hsp26+3xPHO/g" -e "s/DPRE2hsp26_/Hsp26/g" -e "s/DWT_/WT/g")    
    #wc -l $file
    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

exit

# Panel f Control
outFile=scoreMapsk250kexp500_r3000bp_NetA_Fig4fControl.tsv
rm -fvr ${outFile}

files="scoreMapsk250kexp500_r3000bp_Fig5_larvae_DWT_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26Pho_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26_pre2NetA_pre1NetA.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2_pre2NetA_pre1NetA.tsv"

for file in $(ls -1 ${files}) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2NetA//g" -e "s/GAGA//g" -e "s/DPRE2_/DPRE2/g" -e "s/DPRE2hsp26Pho_/Hsp26+3xPHO/g" -e "s/DPRE2hsp26_/Hsp26/g" -e "s/DWT_/WT/g")    

    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

# Panel f Control
outFile=scoreMapsk250kexp500_r3000bp_en_Fig4fControl.tsv
rm -fvr ${outFile}

files="scoreMapsk250kexp500_r3000bp_Fig5_larvae_DWT_pre2en_pre1en.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26Pho_pre2en_pre1en.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2hsp26_pre2en_pre1en.tsv scoreMapsk250kexp500_r3000bp_Fig5_LD_DPRE2_pre2en_pre1en.tsv"

for file in $(ls -1 ${files}) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2en//g" -e "s/GAGA//g" -e "s/DPRE2_/DPRE2/g" -e "s/DPRE2hsp26Pho_/Hsp26+3xPHO/g" -e "s/DPRE2hsp26_/Hsp26/g" -e "s/DWT_/WT/g")    

    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

conda run -n DEseq2 Rscript scripts_clean/one-way_ANOVA_Fig4f.R 0
