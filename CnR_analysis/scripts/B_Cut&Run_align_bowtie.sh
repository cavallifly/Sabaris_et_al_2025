#!/bin/bash
if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing fq.gz files and 2) SAM output folder"
	exit 1
fi

echo "dm6 bowtie2 alignment for Cut&Run will start"
for FILE in $1*".fq.gz"; do
	if [[ "$FILE" = *"_2.fq.gz" ]]; then
		continue
	fi
	fq=${FILE##*/}
	name=${fq/".fq.gz"/}
	sam=$2$name".sam"

	if [[ "$FILE" = *"_1.fq.gz" ]]; then
		sam=${sam/_1.sam/.sam}
		report=${name/_1/}".bowtiereport.txt"
		cmd="bowtie2 -p 40 -x /home/michael.szalay/programs/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 "$FILE" -2 "${FILE/_1.fq.gz/_2.fq.gz}" -S "$sam" 2>"$2$report""
#		err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"${name/_1/}"_bowtie.err"
#		out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"${name/_1/}"_bowtie.out"
	else
		report="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$name".bowtiereport.txt"
		cmd="bowtie2 -p 12 -x /home/michael.szalay/programs/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -U "$FILE" -S "$sam" >> "$report" 2>&1"
#		err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$name"bowtie.err"
 #     		out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/"$name"bowtie.out"
	fi

        if [ -f "$sam" ]; then
                echo "$sam already exists"
                continue
        fi
	echo "$cmd"
	$cmd
#	sbatch -p computepart -c 8 -t 1-00:00:00 --mem=8000 --job-name="GS_bwt2_C&R" -e $err -o $out --wrap="$cmd"
done
