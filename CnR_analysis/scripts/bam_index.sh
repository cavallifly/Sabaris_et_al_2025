#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing sorted bam files 2) output folder for bam index files"
	exit 1
fi

echo "bam sorted to bam indexed will start"
for FILE in $1*"_s.bam"; do
	filename=${FILE##*/}
	filename=${filename/"_s.bam"/}
	output=$2$filename"_s.bai"


	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

        if [[ "$FILE" = *"_s.bai" ]]; then
                echo "$FILE already exists"
                continue
        fi

	cmd="samtools index $FILE $output"
$cmd
        echo "$cmd"
    #    err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_bamindex.err"
     #   out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_bamindex.out"
      #  sbatch -p computepart -c 2 -t 2:00:00 --mem=8000 --job-name="GS_bam_index" -e $err -o $out --wrap="$cmd"
done
