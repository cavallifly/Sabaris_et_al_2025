#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing bam files and 2) output folder for text files"
	exit 1
fi

echo "sambamba dedup"
for FILE in $1*"_s.bam"; do
	filename=${FILE##*/}
	filename=${filename/".bam"/}
	output=$2$filename"_dedup.bam"

	if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

	cmd="sambamba markdup -r --hash-table-size 500000 --overflow-list-size 500000 $FILE $output"
        echo "$cmd"
        #err=$2${filename}"sambamba.err"
        #out=$2${filename}"sambamba.out"
	$cmd
	#sbatch -p computepart -c 4 --mem=3200 -e $err -o $out --wrap="$cmd" --job-name="dedup_sambamba" --time="2:00:00"
done
