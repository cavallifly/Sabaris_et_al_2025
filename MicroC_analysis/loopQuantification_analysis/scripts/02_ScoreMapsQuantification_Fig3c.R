library(misha)
library(shaman)
library(doParallel)

### misha working DB
mDBloc <-  "/zdata/data/mishaDB/trackdb/"
db <- "dm6"
dbDir <- paste0(mDBloc,db,"/")
gdb.init(dbDir)
gdb.reload()

source("./scripts_clean/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

k    <- 250
kexp <- 2*k
pTarget <- "ED_PH18"
k    <- 250
kexp <- 2*k
zmin <- -100
zmax <-  100
scoreData <- list(
        ED_PH18       = paste0("microC.mCScores_ED_PH18_merge_dm6_BeS.mCScores_ED_PH18_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb"),
        ED_PH29       = paste0("microC.mCScores_ED_PH29_merge_dm6_BeS.mCScores_ED_PH29_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb")
        )

scoreTracks <- scoreData
print(scoreTracks)

# Genomic locations of interest
pre1dac  = gintervals("chr2L",16422514,16422515)
pre2dac  = gintervals("chr2L",16485929,16485930)

#pre1en   = gintervals("chr2R",11476251,11476252)
#pre2en   = gintervals("chr2R",11528751,11528752)

pre1NetA = gintervals("chrX",14656251,14656252)
pre2NetA = gintervals("chrX",14751251,14751252)
 
#pre1Fab7 = gintervals("chr3R",16899628,16899629)
#pre2Fab7 = gintervals("chr3R",16900158,16900159)
	    
viewpoints <- list(
	   pre1dac  = gintervals("chr2L",16422514,16422515),
	   pre2dac  = gintervals("chr2L",16485929,16485930),
	   #pre1en   = gintervals("chr2R",11476251,11476252),
	   #pre2en   = gintervals("chr2R",11528751,11528752),
	   #pre1NetA = gintervals("chrX" ,14656251,14656252),
	   #pre2NetA = gintervals("chrX" ,14751251,14751252),
	   #pre1Fab7 = gintervals("chr3R",16899628,16899629),
	   #pre2Fab7 = gintervals("chr3R",16900158,16900159)	   
	   )
vPoints <- do.call(rbind.data.frame, viewpoints)

####################### Quantify Significance
vPairs <- build2D(vPoints,vPoints,diag=FALSE, minDist=10e3)

scoreTracks <- scoreData

res <- 3.0e3
testedIntrvs <- c()
pointData <- list()
for(i in 1:nrow(vPoints))
{
    for(j in 1:nrow(vPoints))
    {

        if(i == j){next;}
        if(vPoints[i,"start"] - vPoints[j,"start"] < 10e3){next;}
	if(vPoints[i,"chrom"] != vPoints[j,"chrom"]){next;}

    	intrv <- gintervals.2d(vPoints[i,"chrom"],vPoints[i,"start"],vPoints[i,"end"],vPoints[j,"chrom"],vPoints[j,"start"],vPoints[j,"end"])
	pointName <- paste(rownames(vPoints)[i],rownames(vPoints)[j],sep="_")
	#pdf(paste0("PRE_loop_quantifications_k",k,"_kexp",kexp,"_Fig4_",pTarget,"_",pointName,".pdf"),width=14)

	#layout(matrix(1:3,ncol=3,byrow=TRUE),respect=FALSE)

	testedIntrvs <- rbind(testedIntrvs,intrv)

	print(paste(rownames(vPoints)[i],rownames(vPoints)[j]))
	print(pointName)
	print(intrv)
	intrv <- expand2D(intrv,res)
	print(intrv)
	for(set in names(scoreTracks))
	{
	    setName = gsub('hic.','',gsub('_merge_dm6_SD','',gsub('_merge_dm6_BeS','',gsub('hicScores_','',set))))
	    print(setName)
	    data <- gextract(scoreTracks[[set]],intrv,colnames="score")
	    pointData[[pointName]][[setName]] <- data$score
	    outMatrix <- paste0("scoreMapsk",k,"kexp",kexp,"_r",res,"bp_Fig4_",set,"_",pointName,".tsv")
	    write.table(data[,c(1,2,3,4,5,6,7)], file = outMatrix, sep="\t", row.names=FALSE, quote=FALSE)
	}
	#jitterPlot(pointData[[pointName]],pTarget=pTarget,main=pointName, ylim=c(0,100),jitWidth=.75,cols=c("darkslateblue",rep("darkslategrey",length(names(scoreTracks))-1)))
	####

	#dev.off()	
    }
}
