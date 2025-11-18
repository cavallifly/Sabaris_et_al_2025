#!/bin/bash

#SBATCH --job-name 06_accepted
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 06_get_accepted_loops.out # File to which STDOUT will be written
#SBATCH -e 06_get_accepted_loops.out # File to which STDERR will be written

assembly=dm6

for inDir in $(ls -1 | grep Loops | grep resD | grep filtered_merged) ;
do
    if [[ ! -d $inDir ]];
    then
	continue
    fi
    echo $inDir

    cd $inDir
    outFile=${inDir}_accepted_merged.bedpe
    
    ls -1 ./acceptedLoops/*png | sed -e "s,/, ,g" | awk '{print $NF}' | sed "s/_/ /g" | awk '{printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$3,$4,$5,$6,$7)}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4n -k 5,5n -k 6,6n | python ~/PyGLtools/merge.py -stdInA -d 1 | grep -v chrA | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4n -k 5,5n -k 6,6n > ../${outFile}

    cd ..

    continue    
    check=$(docker exec higlass-container python higlass-server/manage.py list_tilesets | grep "${outFile%.bedpe}" | wc -l)
    if [[ $check -eq 1 ]];
    then
	echo "${outFile} already loaded in higlass! Going to the next one!"
	continue
    fi

    continue
    # Prepare bedpe files for higlass upload
    conda run -n higlass clodius aggregate bedpe --assembly dm6 \
	--chr1-col 1 --from1-col 2 --to1-col 3 --chr2-col 4 --from2-col 5 --to2-col 6 \
	   --output-file ${outFile%.bedpe}.loops \
	   ${outFile}           
    
    # Upload loops in higlass
    docker exec higlass-container ls /tmp
    rsync -avz ${outFile%.bedpe}.loops /zdata/data/hg-tmp
    docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/${outFile%.bedpe}.loops --filetype bed2ddb --coordSystem dm6 --datatype 2d-rectangle-domains --project-name epiCancer --name "${outFile%.bedpe}"
    rm -fvr /zdata/data/hg-tmp/${outFile%.bedpe}.loops
    
done # Close cycle over $condition
