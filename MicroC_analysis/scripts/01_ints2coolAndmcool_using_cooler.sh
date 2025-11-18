#!/bin/bash

targetResolution=$1
if [[ ${targetResolution} != "" ]]
then
    echo "You choose to write in a .tsv format the matrices at ${targetResolution} bp"
    outTextDir=balanced_text_matrices
    mkdir -p ${outTextDir}
fi

outDir=01_cool_files
mkdir -p ${outDir}
cd ${outDir}

inFiles=$(ls -1 ../01_ints_files/*Rep*ints | grep -v OLD | grep -v DS)
echo $inFiles
resolution=1000

for inFile in $inFiles ;
do
    #(
	echo $inFile
	
	assayName=$(    echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $1}')
	targetName=$(   echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $2}')
	cellName=$(     echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $3}')
	condition=$(    echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $4}')
	replicateName=$(echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $5}')
	laneName=$(     echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $6}')
	assembly=$(     echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $7}')
	tag=$(          echo ${inFile} | sed "s,/, ,g" | awk '{print $NF}' | sed "s,_, ,g" | awk '{print $8}')
	
	chrSizes=/home/common_pipelines/utils/chrom_sizes_${assembly}_higlass.txt
	#cat ${chrSizes}

	outName=$(echo $inFile | sed "s,/, ,g" | awk '{print $NF}' | sed -e "s/\.allValidPairs//g" -e "s/\.ints//g")
	
	resolution=1000
	if [[ ${assayName} == "microc" || ${assayName} == "mC" ]] ;
	then
	    resolution=100 # The .cool file is generated at 100 bp for microc samples
	fi
	if [[ ${assayName} == "hic" || ${assayName} == "chic" ]] ;
	then
	    resolution=1000 # The .cool file is generated at 100 bp for microc samples
	fi        
	
	outFile=${outName}_${resolution}bp.cool
	echo $outFile
	if [[ ! -e ${outFile} ]];
	then
	    touch ${outFile}
	    echo $inFile $outFile	
	    awk '{if($1=="chrom1"){next}; for(i=0;i<$7;i++){print $0}}' $inFile > _tmp_${outFile}
	    conda run -n cooler cooler cload pairs --assembly ${assembly} --chrom1 1 --pos1 2 --chrom2 4 --pos2 5 ${chrSizes}:${resolution} _tmp_${outFile} ${outFile}
	    rm -fvr _tmp_${outFile}
	fi
	if [[ ! -s ${outFile} ]];
	then
	    continue
	fi
	
	# A second .cool file is generated at 20000 bp coarsening the first one for the hicrep analysis
	hicrepRes=20000
	coarsening=$(echo ${hicrepRes} $resolution | awk '{print int($1/$2)}')
	
	inFile=${outFile}
	outFile=${outName}_${hicrepRes}bp.cool
	echo $inFile $outFile
	
	if [[ ! -e ${outFile} ]];
	then
	    touch ${outFile}	
	    conda run -n cooler cooler coarsen -k ${coarsening} ${inFile} -o ${outFile}
	    #-k, --factor INTEGER     Gridding factor. The contact matrix is
            #                         coarsegrained by grouping each chromosomal contact
            #                         block into k-by-k element tiles  [default: 2]
	fi
    #) &
done # Close cycle over $inFile
wait

echo
echo "# Merging the replicates"
for mergeCoolFile in $(ls -1 *_Rep1_*bp.cool | grep -v 20000bp | sed "s,_Rep1_,_merge_,g");		     
do
    echo $mergeCoolFile

    assayName=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $1}')
    targetName=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $2}')
    cellName=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $3}')
    condition=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $4}')
    replicateName=$(echo Rep)
    laneName=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $6}')
    assembly=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $8}')
    resolution=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $9}')
    
    if [[ ${assayName} != "chic" ]];
    then
	resolution=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $9}')

	replicatesName=${assayName}_${targetName}_${cellName}_${condition}_?ep*_${laneName}_${assembly}_${tag}_${resolution}
    fi

    if [[ ${assayName} == "chic" ]];
    then
	onTarget=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $9}')
	resolution=$(echo ${mergeCoolFile} | sed "s,_, ,g" | awk '{print $10}')

	replicatesName=${assayName}_${targetName}_${cellName}_${condition}_?ep*_${laneName}_${assembly}_${tag}_${onTarget}_${resolution}
    fi
    echo ${replicatesName}
    
    nReplicates=$(ls -1 ${replicatesName} | wc -l)
    if [[ ${nReplicates} -lt 2 ]];
    then
	echo "Only one replica avalilable for"
	ls -1 ${replicatesName}
	echo "I am skipping merging"
	echo
	continue
    fi

    coolFiles=$(ls -1 ${assayName}_${targetName}_${cellName}_${condition}_?ep*_all_${assembly}_${tag}_${resolution})
    echo "Merging ${coolFiles} into ${mergeCoolFile}"
    if [[ ! -e ${mergeCoolFile} ]];
    then
	conda run -n cooler cooler merge ${mergeCoolFile} ${coolFiles}
    fi
    echo
    
done # Close cycle over $coolFile
#exit

echo "Creating mcool files"
for coolFile in $(ls -1 *bp.cool | grep merge | grep microc | grep -v 20000bp);
do
    echo $coolFile
    #conda run -n cooler cooler info ${coolFile} &>> coolFiles_info.out
    
    assayName=$(echo ${coolFile} | sed "s,_, ,g" | awk '{print $1}')

    if [[ ${assayName} == "microc" ]];
    then
	#resolutions="100,200,400,600,800,1000,2000,4000,8000,10000,20000,50000,100000,1000000,5000000,10000000" # As in Paldi et al.
	resolutions="100,200,400,800,1000,2000,4000,8000,10000,20000,40000,80000,100000,200000"
    else
	#resolutions="1000,2000,5000,10000,20000,50000,100000,500000,1000000,2000000,5000000,10000000"
	resolutions="1000,2000,4000,8000,10000,20000,40000,80000,100000,200000"
    fi

    mcoolFile=$(echo ${coolFile} | sed -e "s,\.cool,\.mcool,g" -e "s,_1000bp,,g" -e "s,_100bp,,g")
    echo $mcoolFile    
    if [[ ! -e ${mcoolFile} ]];
    then
	touch ${mcoolFile}
    
	echo "cool2mcool: ${coolFile} ${mcoolFile} containing the resolutions ${resolutions}"
	conda run -n cooler cooler zoomify --balance -r ${resolutions} ${coolFile} -o ${mcoolFile} &>> cool2mcool_${coolFile%.cool}.out
    fi
    
    for res in $(echo ${resolutions} | sed "s,\,, ,g");
    do
	echo $mcoolFile >> mcoolFiles_info.out
	#conda run -n cooler cooler info ${mcoolFile}::/resolutions/${res} >> mcoolFiles_info.out 2> /dev/null

	if [[ ${res} -eq ${targetResolution} ]];	   
	then
	    echo $res
	    conda run -n cooler python /home/common_pipelines/microc/03_cooltools_compute_marginals.py ${targetResolution}   

	    outTextFile=${outTextDir}/${mcoolFile%.mcool}_${targetResolution}bp.tsv

	    if [[ ! -e ${outTextFile} ]];
	    then
		awk 'BEGIN{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","chrom1","start1","end1","chrom2","start2","end2","count","ICE")}' > ${outTextFile}   
		cooler dump -H --na-rep "NA" --balanced --join ${mcoolFile}::/resolutions/${res} | grep -v chrom1 | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > ${outTextFile}_tmp
		awk '{if(NF==4){bin=$1"_"$2"_"$3; bad[bin]=1}else{bin1=$1"_"$2"_"$3; bin2=$4"_"$5"_"$6; if(bad[bin1]==1 || bad[bin2]==1){next}else{print $0}}}' ${outTextDir}/MADmax_badColumns_${mcoolFile%.mcool}_at_${res}bp.tsv ${outTextFile}_tmp >> ${outTextFile}
		rm -fvr ${outTextFile}_tmp		
	    fi
	    grep -i na ${outTextFile}
	fi       
    done # Close cycle over resolution
    echo
	
done # Close cycle over coolFile

# Generating the .mcool files


cd ..
