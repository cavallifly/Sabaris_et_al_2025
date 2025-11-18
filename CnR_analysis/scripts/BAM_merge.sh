#!/bin/bash
bam1=/zdata/data/gsabaris/2024_GAF/bam/E_PH_R3.bam
bam2=/zdata/data/gsabaris/2024_GAF/bam/E_PH_R2.bam
bam3=/zdata/data/gsabaris/2024_GAF/bam/E_PH_R1.bam
#bam4=/zdata/data/gsabaris/2024_GAF/bam/E_H3K27me3_R1.bam
output=/zdata/data/gsabaris/2024_GAF/bam/E_PH_R1R2R3_merged.bam
	
	cmd="samtools merge $output $bam1 $bam2 $bam3 $bam4"
        echo "$cmd"
$cmd
#        err="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/merge_bam.err"
 #       out="/work/cavalli/gsabaris/Cut_and_Run/slurmlog/merge_bam.out"
  #      sbatch -p computepart -c 8 -t 1:00:00 --mem=1200M --job-name="GS_merge_bam" -e $err -o $out --wrap="$cmd"
