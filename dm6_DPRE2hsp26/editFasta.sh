mutation=dm6_DPRE2hsp26

### Deletion
#Coordinates of the deletion: chr2L:16485929-16486572
dTargetChrom="chr2L"
dStart=16485929 # Starting nt of the deletion (included)
dEnd=16486572  # Ending   nt of the deletion (included)
dLength=$((${dEnd}-${dStart}+1))
echo "Deletion ${dTargetChrom}:${dStart}-${dEnd} Length ${dLength}"
#Deleted sequence:


###Replacement
echo "Deleted region will be replaced by a sequence of ${dLength} n (undetermined nt)."
echo "This is handy when we want to refer to the insertion, which will have the"
echo "same coordinates as in the WT genome"
#replaceSeq=$(echo $dLength | awk '{for(i=0;i<$1;i++)printf("%s","n")}')
replaceSeq="" #$(echo $dLength | awk '{for(i=0;i<$1;i++)printf("%s","n")}')
echo $replaceSeq
echo

### Insert
#Coordinates of the insertion: chr2L:16560363-
iTargetChrom="chr2L"
iStart=16485929
insertSeq="acccgctggagcttttctagaagagtccggaaagtgacagaaaaaggaaagagtagagagagagaagagaagagagagaacgtgcacagagagaaaaaattatgaatatagctttgtatcgtttgataaacgagggatcaaattagaaatttcttttataaaagtctaaattgaatgctgctctaaaaacaaaaaataaaactaactaacctttcccaataaatgccatgagtttcgagaaattttaaatgtagttgttgtcgtacgcgcaaaagcagagctagtgcgcctgtatgagtgagagagccgaagtttctagagagttgttgcaattttctagaaagagcgcaaaagaaacagccggctaa"
iLength=$(echo ${insertSeq} | awk '{printf("%s",$1)}' | wc | awk '{print $NF}')
replaceSeq=$(echo $dLength ${iLength} | awk '{for(i=$2;i<$1;i++)printf("%s","n")}')
insertSeq=${insertSeq}${replaceSeq}
iLength=$(echo ${insertSeq} | awk '{printf("%s",$1)}' | wc | awk '{print $NF}')
iEnd=$((${iStart}+${iLength}-1))
echo "Insertion ${iTargetChrom}:${iStart}-${iEnd} Length ${iLength}"
#NOTE: the inserted sequence is the PRE1 sequence +/ 50nt, that are the result of the recombination process
echo ${insertSeq}

inFile=/zdata/data/DBs/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.fa
lLength=$(head -2 $inFile | tail -1 | awk '{print(length($0))}')
outFile=${mutation}.fa
rm -fr ${outFile}
echo

#Visualize the changed regions: chr2L:16420000-16570000
echo "-> Difference should be detected in $((${dStart}-16420000))-$((${dStart}-16420000)) for the replacement at PRE1"
echo "-> Difference should be detected in $((${dEnd}-16420000))-$((${dEnd}+${dLength}-16420000)) for the insertion"
vTargetChrom="chr2L"
vStart=16420000
vEnd=16570000
echo "To Blast the INDELs we write the region ${vTargetChrom}:${vStart}-${vEnd} for the original and the modified genome"
echo

for chr in $(grep ">" $inFile | sed "s/>//g")
do
    echo $chr
    if [[ -e ${chr}.seq ]];
    then
	cat ${chr}.seq >> ${outFile}
	continue
    fi
    
    fLine=$(grep -n ${chr} ${inFile} | sed "s/:/ /g" | awk '{print $1}')
    lLine=$(grep -n ">" ${inFile} | sed "s/:/ /g" | awk -v fLine=${fLine} '{if($1>fLine){print $1-1; exit}}')
    if [[ $lLine == "" ]];
    then
	lLine=$(wc -l $inFile | awk '{print $1}')
    fi
    echo $fLine $lLine

    # If the chromosome is not the target file, just write it in the new fasta file
    if [[ $chr != ${dTargetChrom} ]];
    then
	awk -v fLine=${fLine} -v lLine=${lLine} '{if(fLine <= NR && NR <= lLine){print $0}}' ${inFile} >> ${outFile}

	#Write the chromosome sequence
	outFileChr=${chr}.seq
	awk -v fLine=${fLine} -v lLine=${lLine} '{if(fLine <= NR && NR <= lLine){printf("%s\n",$0)}}' ${inFile} > ${outFileChr}
    fi

    if [[ $chr == ${dTargetChrom} ]];
    then
	#Store the entire sequence of the ${dTargetChrom}
	targetChromSeq=$(awk -v fLine=${fLine} -v lLine=${lLine} '{if(fLine < NR && NR <= lLine){printf("%s",$0)}}' ${inFile})
	#Get the chromosome length
	targetChromLen=$(echo ${targetChromSeq} | awk '{printf("%s",$1)}' | wc | awk '{print $NF}')
	echo $targetChromLen

	#Do the deletion
	#Write the sequence of ${dTargetChrom} from start (0) to the nucleotide right before the deletion (${dStart}-1)
	echo ${targetChromSeq} | awk -v dStart=${dStart} '{print substr($0,1,dStart-1)}' > _tmp
	#Write the sequence of ${dTargetChrom} from the nucleotide right after the deletion (${dEnd}+1) to the end (N)
	echo ${targetChromSeq} | awk -v dStart=${dStart} -v dEnd=${dEnd} '{print substr($0,dEnd+1)}' >> _tmp
	newTargetChromSeq=$(awk '{printf("%s",$1)}' _tmp)
	rm _tmp
	
	#Do the insertion
	#Write ${newTargetChromSeq} from start (0) to the nucleotide right before the insertion (${iStart}-1)
	echo ${newTargetChromSeq} | awk -v iStart=${iStart} '{print substr($0,1,iStart-1)}' > _tmp
	echo "Length of the new chromosome before the insertion $(awk '{printf("%s",$1)}' _tmp | wc | awk '{print $NF}')"
	#Write the inserted sequence of the PRE1 +/- 50 nt
        echo ${insertSeq} >> _tmp
	echo "Length of the new chromosome after the insertion  $(awk '{printf("%s",$1)}' _tmp | wc | awk '{print $NF}')"
	#Write ${newTargetChromSeq} from the nucleotide right after the start of the insertion (${iStart}+1) to ${targetChromLen}
	echo ${newTargetChromSeq} | awk -v iStart=${iStart} -v targetChromLen=${targetChromLen} -v iLength=${iLength} '{print substr($0,iStart,targetChromLen-iLength-iStart+1)}' >> _tmp
	#echo ${newTargetChromSeq} | awk -v iStart=${iStart} -v targetChromLen=${targetChromLen} '{print substr($0,iStart)}' >> _tmp	
	newTargetChromLen=$(awk '{printf("%s",$1)}' _tmp | wc | awk '{print $NF}')
	echo "Length of the new chromosome after adding the rest of the chromosome from the insertion point ${newTargetChromLen}"
	echo "Length of the before deletion and replacement ${targetChromLen}"
	echo "To preserve the total length of the chromosome the last $(echo ${targetChromLen} ${newTargetChromLen} | awk '{print $1-$2}') nt have been cut at the telomere"
	#NOTE: We cut the last ${iLength} nucleotides of the chromosome to have the same chromosome sizes
	modTargetChromSeq=$(awk '{printf("%s",$1)}' _tmp)
	rm _tmp

	echo ">${chr}" >> ${outFile}	
	echo ${modTargetChromSeq} | awk -v lLength=${lLength} '{for(i=1;i<=length($0);i+=lLength){printf("%s\n",substr($0,i,lLength))}}' >> ${outFile}

	outFileBlast=${chr}_original.fasta
	echo ">${chr}" > ${outFileBlast}
	echo ${targetChromSeq}    | awk -v vStart=${vStart} -v vEnd=${vEnd} -v lLength=${lLength} '{for(i=vStart;i<=vEnd;i+=lLength){printf("%s\n",substr($0,i,lLength))}}' >> ${outFileBlast}

	outFileBlast=${chr}_modified.fasta
	echo ">${chr}" > ${outFileBlast}
	echo ${modTargetChromSeq} | awk -v vStart=${vStart} -v vEnd=${vEnd} -v lLength=${lLength} '{for(i=vStart;i<=vEnd;i+=lLength){printf("%s\n",substr($0,i,lLength))}}' >> ${outFileBlast}

	#Write the chromosome sequence
	outFileChr=${chr}.seq
	awk -v fLine=${fLine} -v lLine=${lLine} '{if(fLine <= NR && NR <= lLine){printf("%s\n",$0)}}' ${inFile} > ${outFileChr}
	
    fi

done

# Building new bowtie2 index
#bowtie2-build ${mutation}.fa ${mutation}
#exit

# Saving it in the database
mkdir -p /home/minHiC2_sourceFiles/bowtie2/idxs/${mutation}/
rsync -avz ${mutation}.* /home/minHiC2_sourceFiles/bowtie2/idxs/${mutation}
ls -lrtha /home/minHiC2_sourceFiles/bowtie2/idxs/${mutation}/

# Preparing the mishaDB entry
mkdir -p /zdata/data/mishaDB/trackdb/${mutation}/seq/
rsync -avz *seq /zdata/data/mishaDB/trackdb/${mutation}/seq/
sed "s/XXXmutationXXX/${mutation}/g" scripts/chrom.keys > chrom.keys
mkdir -p /zdata/data/mishaDB/trackdb/${mutation}/seq/redb/
mkdir -p /zdata/data/mishaDB/trackdb/${mutation}/tracks/
#exit

rsync -avz /zdata/data/mishaDB/trackdb/dm6/chrom_sizes.txt /zdata/data/mishaDB/trackdb/${mutation}/
sed "s/XXXmutationXXX/${mutation}/g" scripts/mapa_template.config   > mapa_${mutation}.config
sed "s/XXXmutationXXX/${mutation}/g" scripts/redb_template.conf     > redb_${mutation}.conf
sed "s/XXXmutationXXX/${mutation}/g" scripts/create_redb_template.R > scripts/create_redb_${mutation}.R
Rscript scripts/create_redb_${mutation}.R &> create_redb_${mutation}.out
exit
