
for sample in Embryo_GAGA14 Embryo_GAGA34 Embryo_GAGAmut Embryo_WTGAGA ;
do
    Rscript ./scripts/PEscan_parallel.R $sample &>> PEscan_parallel_${sample}.out &
done
