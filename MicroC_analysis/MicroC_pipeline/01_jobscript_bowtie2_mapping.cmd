
### USED  ###
bowtieIndex=DBs/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome

bowtie2 --very-fast --no-unal -k 1 -p 10 -x $bowtieIndex -1 ./mC_${sample}/mC_${sample}_L1_1.fq.gz -2 ./mC_${sample}/mC_${sample}_L1_2.fq.gz -S ./mC_${sample}/mC_${sample}_L1.sam &>> microC_Logs_${sample}.txt

