#!/bin/bash
if [ "$#" -ne 4 ]; then
        echo "Please specify 1) IP rep1 file dir/index 2) IP rep2 file dir/index 3) INPUT file dir/index and 4) OUTPUT folder"
	exit 3
fi
echo "macs3 callpeak will starts"

if [ "$#" -eq 4 ]; then
        for FILE in $1*"dedup.bam"; do
        file1=$1*"dedup.bam"
	file2=$2*"dedup.bam"
	INPUT=$3*"dedup.bam"
        q="0.01"
	output=${FILE##*/}
        output=${output/"_R1_s_dedup.bam"/"_R1R2merged_macs3_q"$q}

        if [ -f "$output" ]; then
                echo "$output already exists"
                continue
        fi

	cmd="macs3 callpeak -t $file1 $file2 -c $INPUT -g dm -n $output -f BAMPE -q $q --outdir $4"
        echo $cmd
	#err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/${output##*/}.err"
        #out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/${output##*/}.out"
	$cmd
        #sbatch -c 12 -p computepart -t 2:00:00 --mem=32000 --job-name="GS_macs2" -e $err -o $out --wrap="$cmd"
	done
fi
