
for sample in Embryo WD PH18 ;
do
    Rscript ./scripts/01_PEscan_parallel.R $sample &>> PEscan_parallel_${sample}.out &
done
