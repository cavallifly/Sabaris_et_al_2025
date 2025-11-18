#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing bam files and 2) output folder for sorted _s.bam files"
	exit 1
fi

echo "sort bam files"
for FILE in $1*".bam"; do
	filename=${FILE##*/}
	filename=${filename/".bam"/}
	output=$2$filename"_s.bam"

	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

        if [[ "$FILE" = *"_s.bam" ]]; then
                echo "$FILE already exists"
                continue
        fi

	cmd="samtools sort $FILE -o $output"
        echo "$cmd"
#        err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_sort_bam.err"
 #       out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_sort_bam.out"
  	$cmd
#      sbatch -p computepart -c 1 -t 2:00:00 --mem=8000 --job-name="GS_sortBAM" -e $err -o $out --wrap="$cmd"
done
