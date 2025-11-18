
# Panel c
outFile=scoreMapsk250kexp500_r3000bp_dac_Fig3c.tsv
rm -fvr ${outFile}
for file in $(ls -1 *PH*dac*tsv) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/PH18_/Control/g" -e "s/PH29_/PH_KD/g" -e "s/pre2dac//g")
    #wc -l $file
    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done

awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

conda run -n DEseq2 Rscript scripts_clean/one-way_ANOVA_Fig3c.R 0
