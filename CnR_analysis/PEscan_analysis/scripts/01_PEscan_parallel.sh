
for sample in ED Embryo WD ;
do
    Rscript ./scripts/PEscan_parallel.R $sample &>> PEscan_parallel_${sample}.out &
done
