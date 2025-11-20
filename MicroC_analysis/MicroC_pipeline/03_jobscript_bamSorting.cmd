
threads=10

samtools sort -n ${inBam} -@ $threads -o ${outBam} &>> 03_bamSorting.out
