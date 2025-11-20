minD=200
shif=0
prevID=0
multiCounter=0
alignedPairs=0
selfLigation=0
validPairs=0

echo "### DUPLICATES REMOVAL : START ###"
inInts=$2
outDeDupFile=$1

cat <(cat ${inInts} 2> /dev/null) <(gunzip -c ${inInts} 2> /dev/null) | awk -v outFile=${outDeDupFile} 'BEGIN{printf("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tobs\n") > outFile; totalInts=0; filteredInts=0;}{totalInts++; intID=$1"_"$2"_"$3"_"$4; if(ints[intID]==1){next}else{filteredInts++; print $0 > outFile; ints[intID]=1}}END{print "totalPairs: "totalInts; print "uniquePairs: "filteredInts}' &>> 05_dupRemoval.log

echo "### DUPLICATES REMOVAL : END ###"
