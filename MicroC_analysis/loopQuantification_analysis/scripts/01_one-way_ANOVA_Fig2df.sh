# Panel d
outFile=scoreMapsk250kexp500_r3000bp_dac_Fig2d.tsv
rm -fvr ${outFile}
for file in $(ls -1 *GA*pre2dac*tsv) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2dac//g" -e "s/GAGA//g" -e "s/14_/\(GA\)n_mut4/g" -e "s/34_/\(GA\)n_mut2/g" -e "s/mut_/\(GA\)n_mut6/g" -e "s/WT_/WT/g")
    #wc -l $file
    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}
done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n  
echo

# Panel f
outFile=scoreMapsk250kexp500_r3000bp_NetANetB_Fig2f.tsv
rm -fvr ${outFile}
for file in $(ls -1 *GA*pre2NetA*tsv) ;
do
    name=$(echo $file | sed "s/_/ /g" | awk '{print $5"_"$6}' | sed -e "s/pre2NetA//g" -e "s/GAGA//g" -e "s/14_/\(GA\)n_mut4/g" -e "s/34_/\(GA\)n_mut2/g" -e "s/mut_/\(GA\)n_mut6/g" -e "s/WT_/WT/g")
    #wc -l $file
    echo $file $name    
    grep -v chrom1 $file | awk -v n=$name '{print n,$NF}' >> ${outFile}

done
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n

conda run -n DEseq2 Rscript ./scripts/one-way_ANOVA_Fig2df.R 0
