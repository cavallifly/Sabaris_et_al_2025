#!/usr/bin/Rscript

require(misha)
require(shaman)
require(reshape2)

args = commandArgs(trailingOnly=TRUE)

sample = as.character(args[1])

require(doParallel)
registerDoParallel(cores=1)

### misha working DB
mDBloc <- '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

### Auxiliary functions
source('/zdata/data/auxFunctions/auxFunctions.R')
source('/zdata/data/auxFunctions/shamanPEscan_mod1.R')

### General R options
options(scipen=9, gmax.data.size=1e8, gmultitasking=TRUE)
chrs <- as.character(gintervals.all()$chrom)

### Load TSS and genes coordinates
#tssCoordinates <- gintervals.load('intervals.ucscCanTSS')
#rownames(tssCoordinates) <- tssCoordinates$geneName
#geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
#rownames(geneCoordinates) <- geneCoordinates$geneName

### Load the peaks: chrom start end k, where k is the cluster index
#peakFiles <- list.files(path = ".", pattern = c(".+bedpe"))
#peakFiles <- list.files(path = ".", pattern = c("Loops_ED_PH18.+bedpe"))
#peakFiles <- list.files(path = ".", pattern = c(".+AllC.+AllC.+4000.+bedpe"))
peakFiles <- list.files(path = ".", pattern = c(".+AllE.+AllE.+4000.+bedpe"))
print(peakFiles)

#peakFiles <- paste0("peaks/",list.files(path = "peaks/", pattern = c("WD_GAF.")))  # "E_GAF.+Q1"
#peakFiles <- paste0("../peaks/",list.files(path = "../peaks/", pattern = c("GAF_WD_Q4_tethering_k")))
#peakFiles <- c("cluster_1_peaks.bed","cluster_2_peaks.bed","cluster_3_peaks.bed","cluster_4_peaks.bed","cluster_5_peaks.bed","cluster_6_peaks.bed")
#print(peakFiles)
#quit()

### Define the pairs of cluster to compare
### PARAMETERS
p1 <- c() # peaks2
p2 <- c() # peaks1
p3 <- c() # minDist
p4 <- c() # maxDist
p5 <- c() # range

for(peak1 in 1:length(peakFiles))
{
    for(peak2 in 1:length(peakFiles))
    {    
        if(peak1 != peak2){next} # for compare only within the same peaks
        #if(peak1 > peak2){next} # to compare all vs all

        ### Proximal 1A
        minDist   <-  3e3   # Lower limit of the band to consider
    	maxDist   <- 10e6   # Upper limit of the band to consider	
    	range     <-  4e3   # Region to consider centered around the peak
	p1 <- append(p1,peakFiles[peak1])
	p2 <- append(p2,peakFiles[peak2])
	p3 <- append(p3,minDist)
    	p4 <- append(p4,maxDist)
    	p5 <- append(p5,range)
	next
    }
}
pairs <- cbind(p1,p2,p3,p4,p5)
print(paste0(pairs))

### List of runs to analyze. Each entry of this list is a list including the obs and expected (shuffled) tracks

if(sample == "PH18")
{
    allRuns <- list(mC_ED_PH18 = list("microC.mC_ED_PH18_Rep1_dm6_BeS.mC_ED_PH18_Rep1_dm6_BeS", "microC.mC_ED_PH18_Rep1res_dm6_BeS.mC_ED_PH18_Rep1res_dm6_BeS", "microC.mC_ED_PH18_Rep2_dm6_BeS.mC_ED_PH18_Rep2_dm6_BeS", "microC.mC_ED_PH18_Rep2res_dm6_BeS.mC_ED_PH18_Rep2res_dm6_BeS", "microC.mC_ED_PH18_Rep1_dm6_BeS.mC_ED_PH18_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH18_Rep1res_dm6_BeS.mC_ED_PH18_Rep1res_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH18_Rep2_dm6_BeS.mC_ED_PH18_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH18_Rep2res_dm6_BeS.mC_ED_PH18_Rep2res_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "PH29")
{
    allRuns <- list(mC_ED_PH29 = list("microC.mC_ED_PH29_Rep1_dm6_BeS.mC_ED_PH29_Rep1_dm6_BeS", "microC.mC_ED_PH29_Rep1res_dm6_BeS.mC_ED_PH29_Rep1res_dm6_BeS", "microC.mC_ED_PH29_Rep2_dm6_BeS.mC_ED_PH29_Rep2_dm6_BeS", "microC.mC_ED_PH29_Rep2res_dm6_BeS.mC_ED_PH29_Rep2res_dm6_BeS", "microC.mC_ED_PH29_Rep1_dm6_BeS.mC_ED_PH29_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH29_Rep1res_dm6_BeS.mC_ED_PH29_Rep1res_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH29_Rep2_dm6_BeS.mC_ED_PH29_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PH29_Rep2res_dm6_BeS.mC_ED_PH29_Rep2res_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "PHD11")
{
    allRuns <- list(mC_ED_PHD11 = list("microC.mC_ED_PHD11_Rep1_dm6_BeS.mC_ED_PHD11_Rep1_dm6_BeS", "microC.mC_ED_PHD11_Rep1res_dm6_BeS.mC_ED_PHD11_Rep1res_dm6_BeS", "microC.mC_ED_PHD11_Rep2_dm6_BeS.mC_ED_PHD11_Rep2_dm6_BeS", "microC.mC_ED_PHD11_Rep2res_dm6_BeS.mC_ED_PHD11_Rep2res_dm6_BeS", "microC.mC_ED_PHD11_Rep1_dm6_BeS.mC_ED_PHD11_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PHD11_Rep1res_dm6_BeS.mC_ED_PHD11_Rep1res_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PHD11_Rep2_dm6_BeS.mC_ED_PHD11_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_ED_PHD11_Rep2res_dm6_BeS.mC_ED_PHD11_Rep2res_dm6_BeS_shuffle_500Small_1000High"))
}

if(sample == "PH18DStoPH29")
{
    allRuns <- list(mC_ED_PH18DStoPH29 = list("microC.microc_contacts_ED_PH18DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PH18DStoPH29_merge_all_dm6_BeS", "microC.microc_contacts_ED_PH18DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PH18DStoPH29_merge_all_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "PH29DStoPH29")
{
    allRuns <- list(microc_contacts_ED_PH29DStoPH29 = list("microC.microc_contacts_ED_PH29DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PH29DStoPH29_merge_all_dm6_BeS", "microC.microc_contacts_ED_PH29DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PH29DStoPH29_merge_all_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "PHD11DStoPH29")
{
    allRuns <- list(microc_contacts_ED_PHD11DStoPH29 = list("microC.microc_contacts_ED_PHD11DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PHD11DStoPH29_merge_all_dm6_BeS", "microC.microc_contacts_ED_PHD11DStoPH29_merge_all_dm6_BeS.microc_contacts_ED_PHD11DStoPH29_merge_all_dm6_BeS_shuffle_500Small_1000High"))
}

if(sample == "WD")
{
    # mC_WD
    allRuns <- list(mC_WD = list("microC.mC_WD_WT_Rep1_dm6_BeS.mC_WD_WT_Rep1_dm6_BeS", "microC.mC_WD_WT_Rep2_dm6_BeS.mC_WD_WT_Rep2_dm6_BeS", "microC.mC_WD_WT_Rep3_dm6_BeS.mC_WD_WT_Rep3_dm6_BeS", "microC.mC_WD_WT_Rep1_dm6_BeS.mC_WD_WT_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_WD_WT_Rep2_dm6_BeS.mC_WD_WT_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_WD_WT_Rep3_dm6_BeS.mC_WD_WT_Rep3_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "LD")
{
    # mC_LD
    allRuns <- list(mC_LD = list("microC.mC_LD_WT_Rep1_dm6_BeS.mC_LD_WT_Rep1_dm6_BeS", "microC.mC_LD_WT_Rep2_dm6_BeS.mC_LD_WT_Rep2_dm6_BeS", "microC.mC_LD_WT_Rep3_dm6_BeS.mC_LD_WT_Rep3_dm6_BeS", "microC.mC_LD_WT_Rep1_dm6_BeS.mC_LD_WT_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_LD_WT_Rep2_dm6_BeS.mC_LD_WT_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_LD_WT_Rep3_dm6_BeS.mC_LD_WT_Rep3_dm6_BeS_shuffle_500Small_1000High"))
}
if(sample == "Embryo")
{
    # mC_Embryo
    allRuns <- list(mC_Embryo = list("microC.mC_Embryo_WT_Rep1_dm6_BeS.mC_Embryo_WT_Rep1_dm6_BeS", "microC.mC_Embryo_WT_Rep2_dm6_BeS.mC_Embryo_WT_Rep2_dm6_BeS", "microC.mC_Embryo_WT_Rep3_dm6_BeS.mC_Embryo_WT_Rep3_dm6_BeS", "microC.mC_Embryo_WT_Rep1_dm6_BeS.mC_Embryo_WT_Rep1_dm6_BeS_shuffle_500Small_1000High", "microC.mC_Embryo_WT_Rep2_dm6_BeS.mC_Embryo_WT_Rep2_dm6_BeS_shuffle_500Small_1000High", "microC.mC_Embryo_WT_Rep3_dm6_BeS.mC_Embryo_WT_Rep3_dm6_BeS_shuffle_500Small_1000High"))
}

### Analyze each cluster pairs in parallel
foreach(row=1:nrow(pairs)) %dopar% {
#foreach(row=1:1) %dopar% {
    print(row)
    peak1   <- pairs[row,1]
    peak2   <- pairs[row,2]    
    minDist <- as.numeric(pairs[row,3])
    maxDist <- as.numeric(pairs[row,4])
    range   <- as.numeric(pairs[row,5])
    res     <- 250

    print(paste(peak1,peak2,minDist,maxDist,range,res))

    for(peakFile in c(peak1))
    {
        setName1 <- gsub('../peaks/','',gsub('.bedpe','',peak1))
        setName2 <- gsub('../peaks/','',gsub('.bedpe','',peak2))
	setName <- paste(setName1,setName2,sep='_')
	print(paste0("setname ",setName))
	#	quit()

	peakData <- read.table(peakFile, header=F)
	colnames(peakData) <- c('chrom','start','end')
	print(paste0(setName1,'_peScan.df'))
	save(peakData,file=paste0(setName1,'_peScan.df'))
	head(peakData)

	
	peak2Data <- read.table(peak2, header=F)
	colnames(peak2Data) <- c('chrom','start','end')
	head(peak2Data)

	sourcePeaks <- peakData
        targetPeaks <- peak2Data
	
	ints2d <-  build2D(center1D(sourcePeaks),center1D(targetPeaks),diag=FALSE,retNames=FALSE,minDist=minDist,maxDist=maxDist)
	#peIterator <- ints2d[,1:6]
	peIterator <- read.table(peakFile, header=F)
	colnames(peIterator) <- c('chrom1','start1','end1','chrom2','start2','end2')
	print(peIterator)
        print(nrow(peIterator))

        for(run in names(allRuns))
        {		
            runName <- paste('distPEscan_',run,'_minDist',minDist/1000,'_maxDist',maxDist/1000,sep='')
	    print(paste0("RunName"," ",runName))
	    set <- paste(setName,sep='_')
            if(length(allRuns[[run]]) > 2)
            {
    	        oIdxs <- 1:(length(allRuns[[run]])/2)
	        eIdxs <- ((length(allRuns[[run]])/2)+1):(length(allRuns[[run]]))
	        set1 <- unlist(allRuns[[run]][oIdxs])
	        set2 <- unlist(allRuns[[run]][eIdxs])
	    }else{
	        set1 <- allRuns[[run]][[1]]
	        set2 <- allRuns[[run]][[2]]
	    }
	    print(set1)
	    print(set2)

	    #######################################################
	    outFolder <- '/zdata/data/mdistefano/2022_06_08_Project_on_PREs_contacts/04_loopCalling_analysis/peScan_analysis/'
	    #outFolder <- ''
	    outFolder  <- paste0(outFolder,runName,'/')
	    dataFolder <- paste0(outFolder,'RData/')


	    if(!dir.exists(outFolder )){dir.create(outFolder, mode='7777');}
	    if(!dir.exists(dataFolder)){dir.create(dataFolder, mode='7777');}
		
	    outName <- paste0(runName)
	    print(paste("outName",outName,sep=' '))

            # List of XXX
	    outListName <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_peScan.list')
	    outListNameText <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_peScan_text.list')	    
	    outPairListNameText <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_bedpe_text.list')
	    #next
	    if(file.exists(outPairListNameText))
	    {
	        print(paste("File",outListName,"exists",sep=' '))
		next
	    }
	    write.table(peIterator,file=outPairListNameText,sep="\t",quote=F,row.names=F,append=F,col.names=F)
	    print(paste("outListName",outListName,sep=' '))
	    #quit()
	    		
	    allIntervals <- list()
	    allIntervals[[outName]] <- peIterator
		
	    resPEdata <- list()
	    print(paste("allIntervals",names(allIntervals),sep=' '))
    	    for(j in names(allIntervals))
	    {
		print('RUN')
	    	obs <- set1
	    	exp <- set2

		outFig      <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res,'_peScan.png')
		print(paste("outFig",outFig,sep=' '))
		# Matrix of plotted data at resolution=res
		matName     <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res,'_peScan.mtx')
		matNameText <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res,'_peScan_text.mtx')
		print(paste("matName",outListName,sep=' '))


	    	grids <- list()
	    	for(i in 1:length(obs))
	    	{
		    grids[[i]] <- shaman_generate_feature_grid_2d(peIterator, obs[i], exp[i], range = range, resolution = res, min_dist = minDist, max_dist = maxDist)
	    	}

		# Plot at resolution=resch
	    	resPEdata[[set]][[j]] <- my_shaman_plot_feature_grid(grids, range, res, res,  type = "enrichment", fig_fn = outFig, zlim=c(-1,1), text=TRUE, N=nrow(peIterator))
		save(resPEdata, file=outListName)
		write.table(resPEdata, file=outListNameText)
	    }
		
	    for(cl in names(allIntervals))
	    {
		for(set in names(resPEdata))
	    	{
		    y <- resPEdata[[set]][[cl]]
	            mat <- log2(((y$obs+1)/(sum(y$obs,na.rm=TRUE)+1))/((y$exp+1)/(sum(y$exp,na.rm=TRUE)+1)))
	            save(mat,file=matName)
		    write.table(mat,file=matNameText)
	    	}
	    }

	    # Plot at resolution=2*res
	    #outFig      <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res*2,'_peScan.png')
	    # Matrix of plotted data at resolution=2*res
	    #matName     <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res*2,'_peScan.mtx')
	    #matNameText <- paste0(dataFolder,outName,'_',setName,'_min',minDist/1e3,'_max',maxDist/1e3,'_range',range/1e3,'_N',nrow(peIterator),'_res',res*2,'_peScan_text.mtx')

	    #print('READY TO PRINT LOW RES')
	    #grids <- get(load(outListName))
	    #y <- my_shaman_plot_feature_grid(grids[[1]], range, res, res*2,  type = "enrichment", fig_fn = outFig, zlim=c(-1,1), text=TRUE, N=nrow(peIterator))
	    #mat <- log2(((y$obs+1)/(sum(y$obs,na.rm=TRUE)+1))/((y$exp+1)/(sum(y$exp,na.rm=TRUE)+1)))
	    #save(mat, file=matName)
	    #write.table(mat,file=matNameText)

	}
    }
    return(0)
}
