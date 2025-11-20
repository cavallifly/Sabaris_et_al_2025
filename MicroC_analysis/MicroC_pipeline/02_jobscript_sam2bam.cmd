sample=$1

# Filter reads mapped with quality < 30 (mapq)
samtools view -S -@ 10 -q 30 -F 4 -b -o ${outBam} ${inSam} &>> 

#rm ./mC_${sample}/mC_${sample}_L1.sam
