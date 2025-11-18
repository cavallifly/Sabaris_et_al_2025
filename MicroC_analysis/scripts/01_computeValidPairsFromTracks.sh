

#inFile=01_computeValidPairsFromTracks_genomeWide.out
#if [[ ! -e ${inFile} ]];
#then
#    Rscript ./scripts/01_computeValidPairsFromTracks.R > ${inFile}
#fi
#outFile=count_validPairs_genomeWide.out
#grep -a hic_ ${inFile} | grep -av "chr\|dac" | awk '{if(NF==3) print $0}' | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}
#echo
#outFile=count_validPairs_chr2L.out
#grep hic_ ${inFile} | grep "chr2L chr2L" | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}
#echo
#outFile=count_validPairs_dac.out
#grep hic_ ${inFile} | grep dac           | awk '{printf("%s\t%s\t%s\t%s\n",$2,$3,$4,$5)}' | sed "s/\"//g" > ${outFile}

outFile=obsContacts_per_sample.tab
cat count_validPairs_genomeWide.out | sed "s/_Rep/ /g" | awk '{print $1,$3}' | awk '{s[$1]+=$2}END{for(i in s){printf("%s\t%d\n",i,s[i])}}' > ${outFile}
cat count_validPairs_chr2L.out | sed -e "s/_Rep/ /g" -e "s/hic_/ /g" | awk '{print $1,$NF}' | awk '{h[$1]+=$2}END{for(i in h) print i,h[i]}' > obsContacts_per_sample_chr2L.tab
cat count_validPairs_dac.out | sed "s/_Rep/ /g" | awk '{print $1,$NF}' | awk '{h[$1]+=$2}END{for(i in h) print i,h[i]}' > obsContacts_per_sample_dac.tab 
