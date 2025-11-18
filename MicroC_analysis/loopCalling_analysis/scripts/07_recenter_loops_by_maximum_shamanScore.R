require(reshape2) # cat tabulate to matrix
library(misha)
library(shaman)
library(dplyr)
library(ggplot2)
#library(ggpubr)
#library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
#source('/zdata/data/auxFunctions/auxFunctions.R')
options("scipen"=999)

mDBloc <- '/work/user/mdistefano/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

obsData <- list(
	larvae_DWT = c("hic.hicScores_larvae_DWT_merge_dm6_BeS.hicScores_larvae_DWT_merge_dm6_BeS_score_k250_kexp500_step1Mb"),
        ED_PH18    = c("microC.mCScores_ED_PH18_merge_dm6_BeS.mCScores_ED_PH18_merge_dm6_BeS_score_k250_kexp500_step1Mb"),
        ED_PH29    = c("microC.mCScores_ED_PH29_merge_dm6_BeS.mCScores_ED_PH29_merge_dm6_BeS_score_k250_kexp500_step1Mb"),
        ED_PHD11   = c("microC.mCScores_ED_PHD11_merge_dm6_BeS.mCScores_ED_PHD11_merge_dm6_BeS_score_k250_kexp500_step1Mb"),
        WD_WT      = c("microC.mCScores_WD_WT_merge_dm6_BeS.mCScores_WD_WT_merge_dm6_BeS_chr_all_k250_kexp500_step500kb_skin250kb"),	
        LD_WT      = c("microC.mCScores_LD_WT_merge_dm6_BeS.mCScores_LD_WT_merge_dm6_BeS_chr_all_k250_kexp500_step500kb_skin250kb"),
        Embryo_WT  = c("microC.mCScores_Embryo_WT_merge_dm6_BeS.mCScores_Embryo_WT_merge_dm6_BeS_chr_all_k250_kexp500_step500kb_skin250kb")
	       )
print(obsData)

userSample <- args[1]

for(sample in names(obsData)[grep(userSample,obsData)])
{
    print(sample)
    intervalsFile = list.files("./",pattern=paste0(".*",sample,".*accepted_merged.bedpe"))
    print(intervalsFile)

    intervals     <- read.table(intervalsFile, header=F)
    intervals     <- intervals[,c(1,2,3,4,5,6)]
    colnames(intervals) <- c("chrom1","start1","end1","chrom2","start2","end2")
    print(head(intervals))

    outfile <- gsub("accepted_merged.bedpe","accepted_merged_recentered.bedpe",intervalsFile)
    print(outfile)

    if(!file.exists(outfile))
    {
        #header <- paste(c("chrom1","start1","end1","chrom2","start2","end2","score"),collapse="\t")
    	#write(header,file=outfile,sep="\t",append=T)

        # Get tracks and create virtual tracks: this depends only on the specific sample
        track2D <- obsData[[sample]]
        print(paste0("Analysing ",sample))
	for(r in 1:nrow(intervals))
	{
	    interval <- intervals[r,]
	    interval <- gintervals.2d(interval$chrom1,interval$start1,interval$end1,interval$chrom2,interval$start2,interval$end2)
	    print(interval)
	    
	    #chrIntrv2D <- build2D(tmpIntervals1,tmpIntervals2,names1=tmpIntervals1$interval,names2=tmpIntervals2$interval)
	    data     <- gextract(track2D,intervals=interval,colnames="score")
	    maxData <- max(data$score)	    
	    print(data)
	    print(maxData)

	    data <- data[data$score == maxData,]
	    data <- data              %>%
	        reframe(chrom1 = unique(chrom1), start1=min(start1), end1=max(end1), chrom2 = unique(chrom2), start2=min(start2), end2=max(end2), score = unique(score))
	    print(data)


	    centre1 <- as.integer((data$end1 + data$start1)*0.5)
	    data$start1 <- centre1 - 2000
	    data$end1 <- centre1 + 2000

	    centre2 <- as.integer((data$end2 + data$start2)*0.5)
	    data$start2 <- centre2 - 2000
	    data$end2 <- centre2 + 2000

            write.table(data,file=outfile,sep="\t",quote=F,row.names=F,col.names=F,append=T)
	}
    } # Close cycle over sample
}
