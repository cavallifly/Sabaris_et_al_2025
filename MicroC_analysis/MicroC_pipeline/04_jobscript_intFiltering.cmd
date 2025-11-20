

inBam=$1
outInts=$2
minD=200
shif=0

awk 'BEGIN{printf("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tobs\n")}' > ${outFile}

samtools view ${inBam} | awk -v shift=${shif} -v minD=${minD} -v outFile=${outInts} 'BEGIN{prevID=0; multiCounter=0; alignedPairs=0; selfLigation=0; validPairs=0;}{id=$1; flag=$2; c1=$3; p1=$4; mapq=$5; cigar=$6; c2=$7; if(c2=="="){c2=c1}; p2=$8; dist=$9; if(id==prevID){multiCounter++; prevID=id; next}; alignedPairs++; if(c2==c1 && sqrt((p2-p1)*(p2-p1))<minD){selfLigation++; prevID=id; next}; s1=p1+shift; s2=p2+shift; e1=s1+1; e2=s2+1; if(s2<s1){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",c1,s1,e1,c2,s2,e2,1) > outFile}else{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",c2,s2,e2,c1,s1,e1,1) > outFile}; validPairs++; prevID=id}END{print "doubleIDs: "multiCounter; print "proposedAlignments: "alignedPairs; print "selfLigation: "selfLigation; print "validPairs: "validPairs; }' >> 04_intFiltering.log
