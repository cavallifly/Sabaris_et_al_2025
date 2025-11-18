#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing sam files and 2) output folder for mapq-sorted .bam files"
	exit 1
fi

echo "remove mapq <30"
for FILE in $1*".sam"; do
	filename=${FILE##*/}
	filename=${filename/".sam"/}
	output=$2$filename".bam"

	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

	cmd="samtools view -b -q 30 $FILE -o $output"
        echo "$cmd"
        #err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_mapq_sort.err"
        #out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_mapq_sort.out"
        $cmd
	#sbatch -p computepart -c 1 -t 2:00:00 --mem=8000 --job-name="GS_mapq_sort" -e $err -o $out --wrap="$cmd"
done
