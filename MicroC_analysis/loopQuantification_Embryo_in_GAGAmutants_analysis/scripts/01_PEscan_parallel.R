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
peakFiles <- list.files(path = ".", pattern = c("Loops_Embryo.+bedpe"))
print(peakFiles)
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

if(sample == "Embryo_GAGA14")
{
    allRuns <- list(GAnmut4 = list("microC.mC_Embryo_GAGA14_Rep1_dm6_GS.mC_Embryo_GAGA14_Rep1_dm6_GS", "microC.mC_Embryo_GAGA14_Rep2_dm6_GS.mC_Embryo_GAGA14_Rep2_dm6_GS", "microC.mC_Embryo_GAGA14_Rep1_dm6_GS.mC_Embryo_GAGA14_Rep1_dm6_GS_shuffle_500Small_1000High", "microC.mC_Embryo_GAGA14_Rep2_dm6_GS.mC_Embryo_GAGA14_Rep2_dm6_GS_shuffle_500Small_1000High"))
}

if(sample == "Embryo_GAGA34")
{
    allRuns <- list(GAnmut2 = list("microC.mC_Embryo_GAGA34_Rep1_dm6_GS.mC_Embryo_GAGA34_Rep1_dm6_GS", "microC.mC_Embryo_GAGA34_Rep2_dm6_GS.mC_Embryo_GAGA34_Rep2_dm6_GS", "microC.mC_Embryo_GAGA34_Rep1_dm6_GS.mC_Embryo_GAGA34_Rep1_dm6_GS_shuffle_500Small_1000High", "microC.mC_Embryo_GAGA34_Rep2_dm6_GS.mC_Embryo_GAGA34_Rep2_dm6_GS_shuffle_500Small_1000High"))

}
if(sample == "Embryo_GAGAmut")
{
    allRuns <- list(GAnmut6 = list("microC.mC_Embryo_GAGAmut_Rep1_dm6_GS.mC_Embryo_GAGAmut_Rep1_dm6_GS", "microC.mC_Embryo_GAGAmut_Rep2_dm6_GS.mC_Embryo_GAGAmut_Rep2_dm6_GS", "microC.mC_Embryo_GAGAmut_Rep1_dm6_GS.mC_Embryo_GAGAmut_Rep1_dm6_GS_shuffle_500Small_1000High", "microC.mC_Embryo_GAGAmut_Rep2_dm6_GS.mC_Embryo_GAGAmut_Rep2_dm6_GS_shuffle_500Small_1000High"))
}   
if(sample == "Embryo_WTGAGA")
{
    allRuns <- list(WT = list("microC.mC_Embryo_WTGAGA_Rep1_dm6_GS.mC_Embryo_WTGAGA_Rep1_dm6_GS", "microC.mC_Embryo_WTGAGA_Rep2_dm6_GS.mC_Embryo_WTGAGA_Rep2_dm6_GS", "microC.mC_Embryo_WTGAGA_Rep1_dm6_GS.mC_Embryo_WTGAGA_Rep1_dm6_GS_shuffle_500Small_1000High", "microC.mC_Embryo_WTGAGA_Rep2_dm6_GS.mC_Embryo_WTGAGA_Rep2_dm6_GS_shuffle_500Small_1000High"))
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
	    outFolder <- '/zdata/data/mdistefano/2022_06_08_Project_on_PREs_contacts/04_loopCalling_analysis/peScan_analysis_GAGAmutants/'
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
