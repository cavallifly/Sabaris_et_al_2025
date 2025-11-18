#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing indexeded bam files 2) output folder for bigwig files"
	exit 1
fi

echo "bam to bigWig will start"
for FILE in $1*"_s_dedup.bam"; do
	filename=${FILE##*/}
	filename=${filename/"_s_dedup.bam"/}
	output=$2$filename"_s_dedup.bigWig"


	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

        if [[ "$FILE" = *"_dedup.bigWig" ]]; then
                echo "$FILE already exists"
                continue
        fi

	cmd="bamCoverage --normalizeUsing RPKM --ignoreDuplicates -e 0 -bs 10 -b $FILE -o $output"
	#cmd="bamCoverage --normalizeUsing BPM --ignoreDuplicates -e 0 -bs 10 -b $FILE -o $output"
        echo "$cmd"
        #err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_dedup_bamCoverage.err"
        #out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$filename"_dedup_bamCoverage.out"
 	$cmd       
#sbatch -p computepart -c 8 -t 2:00:00 --mem=8000 --job-name="GS_BAMtoBigWig" -e $err -o $out --wrap="$cmd"
done
