#inFiles="Loops_Embryo_WT_BeS.bedpe Loops_Embryo_WT_Batut2020.bedpe Loops_Embryo_WT_Levo2022.bedpe Loops_Embryo_WT_Pollex2024.bedpe Loops_Embryo_WT_Dolsten2025.bedpe"
#S1: Embryo_BeS
#S2: Embryo_Batut2020
#S3: Embryo_Levo2022
#S4: Embryo_Pollex2024
#S5: Embryo_Dolsten2025

for res in 4000;
do
    inFiles=$(ls -1 Loops_*_renamedLoops_ext_${res}bp.bedpe | grep "Embryo_")
    echo $res $inFiles
    for file in ${inFiles} ;
    do
	head $file
	condition=$(echo $file | sed "s/_/ /g" | awk '{print $4}')
	cat ${file} | awk -v c=$condition '{if(NF<8){print $0,c}else{print $0}}' > _tmp
	mv _tmp ${file}
	head $file
    done
    
    cat ${inFiles} | sed -e "s/loop//g" -e "s/_/ /g" | awk '{print $7}' | sort | uniq > categories_at_${res}bp.txt
    head categories_at_${res}bp.txt

    paste <(cat categories_at_${res}bp.txt) <(sed -e "s/S1/ BeS /g" -e "s/S2/ Batut2020 /g" -e "s/S3/ Levo2022 /g" -e "s/S4/ Pollex2024 /g" -e "s/S5/ Dolsten2025 /g" categories_at_${res}bp.txt | awk '{for(i=1;i<=NF;i++){h[$i]++; if(h[$i]==1){printf("%s-",$i)}}; printf("\n") ; for(i=1;i<=NF;i++){h[$i]=0}}') | awk '{print $1,substr($2,1,length($2)-1)}' | sort -k 4,4 > renamed_categories_at_${res}bp.txt
    head renamed_categories_at_${res}bp.txt
    
    rm -fr Loops_AllE_AllE_renamedLoops_ext_${res}bp_*.bedpe *_tmp.bed
    
    for file in ${inFiles} ;
    do

	cat renamed_categories_at_${res}bp.txt <(cat $file | sed -e "s/_/ /g" -e "s/loop//g" | awk '{print $1,$2,$3,$4,$5,$6,$7,$9}') | awk '{if(NF==2){n[$1]=$2}else{print $1,$2,$3,$4,$5,$6,n[$7],$8}}' > ${file%.bed}_tmp.bed
	head ${file%.bed}_tmp.bed
	
	for category in $(awk '{print $2}' renamed_categories_at_${res}bp.txt);
	do
	    
	    awk -v c=${category} '{if($7==c){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}}' ${file%.bed}_tmp.bed >> Loops_AllE_AllE_renamedLoops_ext_${res}bp_${category}.bedpe
	done
    done    
    
done

for file in $(ls -1 Loops_AllE_AllE_renamedLoops_ext_*bp_*.bedpe);
do
    head $file
    #cat $file | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n | uniq > _tmp ; mv _tmp ${file}
    cat $file | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > _tmp ; mv _tmp ${file}    
    echo
    head $file
done
rm -fr *_tmp.bed
