mkdir -p loopFiles
rsync -avz ../mustache_analysis/

fdr=0.200
for dir in $(ls -1 | grep NPC_WT_ | grep default);
do
    if [[ ! -d $dir ]];
    then
	continue
    fi
    
    cd $dir

    for nclusters in $(seq 1 1 5);
    do
	grep -v chrom1 Loop_chr18_52775574_57266687_chr18_52775574_57266687_at_allResolutions_fdr2000_0.000_fdr5000_${fdr}_fdr10000_${fdr}_fdr15000_${fdr}_fdr20000_${fdr}.tsv > _tmp
        condition=$(echo $dir | sed "s/_/ /g" | awk '{print $2"_"$3"_"}')
	echo $condition

	# Cluster based on anchor 1
	Rscript ../scripts/seqPlots_loops_analysis.R _tmp $condition $nclusters 1
	# Cluster based on anchor 2
	#Rscript scripts/seqPlots_loops_analysis.R Loop_chr18_52775574_57266687_chr18_52775574_57266687_at_allResolutions_fdr2000_0.000_fdr5000_0.120_fdr10000_0.120_fdr15000_0.120_fdr20000_0.120.tsv $condition $nclusters 2    
	# Cluster based on both anchors
	#XXX
    done # Close cycle over $nclusters
    cd ..
done # Close cycle over $dir
