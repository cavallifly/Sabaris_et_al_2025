library("fpc")
library("plyr")
#library("tidyverse")
#library("dplyr")
library("ggplot2")
library("reshape2")
library("seqplots")
library(BSgenome)
options("scipen"=999, max.print=999999)
require(doParallel)
registerDoParallel(cores=2)	

source('/zdata/data/auxFunctions/auxFunctions.R')
args = commandArgs(trailingOnly=TRUE)

chromSizes <- read.table("/home/Programs/chrom_sizes_dm6_higlass.txt")
mainDir <- "/zdata/data/mdistefano/2022_06_08_Project_on_PREs_contacts/04_loopCalling_analysis/seqplot_analysis/"


### Get the list of bigWig files ###
bwDir <- paste0(mainDir,"/bigWigs/")
assembly <- 'dm6'
cellline <- args[1]

ClusteringAnchor <- 0 #

safeRowMeans <- function(x) {
  if (nrow(x) < 2) {
    return(0)   # or return(rowMeans(x)) if 1 row is acceptable
  }
  rowMeans(x,na.rm = T)
}

### Get the list of bigWig files ###
testBW <- c(list.files(bwDir, pattern = 'bw'),list.files(bwDir, pattern = 'bigWig')) # tutti i bigWig
#print(testBW)
#testBW <- testBW[ grep("merge", testBW) ]
print(testBW)
#testBW <- testBW[ grep(cellline, testBW) ]
testBW <- testBW[ grep(paste0(cellline,"_"), testBW) ]
#print(testBW)
#targets <- paste("IgG",sep="|")
#print(targets)
#testBW <- testBW[grep(pattern=targets, testBW)] # Select by target and then by cell line
#testBW <- testBW[-grep("INPUT|GRH|PC|PSC|H3K27ac_WD|H3K27me3_WD", testBW)]
print(testBW)
### Get the list of bigWig files ###
#quit()

### Set the parameters of the heatmap-plots and the clustering ###
#range <- 25e3 
#range <- 2e3 # For ED_PH18 to identify PcG loops
range <- 5e3  # For Embryo_WT and WD_WT and ED_PH18
binSize <- range/100
raster <- TRUE
#meanBins <- 25 # For ED_PH18
meanBins <- 100 # For Embryo_WT and WD_WT
#kMax <- 15 # args[1] # Number of clusters
kMax <- 2 # args[1] # Number of clusters
cols    <- list(peaks=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255))
#cols     <- list(peaks=colorRampPalette(c('grey90','grey15'))(200))
### Set the parameters of the heatmap-plots and the clustering ###


#peakFiles <- list.files(paste0(mainDir,"/loopFiles/"), pattern= paste0(cellline,".*reso.*bedpe"))
peakFiles <- list.files(paste0(mainDir,"/loopFiles/"), pattern= paste0(cellline,".bedpe"))
#peakFiles <- list.files(paste0(mainDir,"/loopFiles/"), pattern= paste0(cellline,".*bedpe"))
#print(peakFiles)
#peakFiles <- peakFiles[grep(paste0(cellline,"_"), peakFiles)]
print(peakFiles)
#quit()

for(condition in c(cellline))
{

    wDir <- paste0("./",condition,"/")
    if(!dir.exists(wDir))
    {
        dir.create(wDir, mode="7777", recursive=TRUE)
    }
    setwd(wDir)

    ### Load the loops ###
    ### Load 1st anchors ###
    for(peakFile in peakFiles)
    {
        #wDir <- gsub(".bedpe","",paste0(basename(peakFile),"_notOptimized"))
	wDir <- gsub(".bedpe","",paste0(basename(peakFile)))	
	print(wDir)
    	if(dir.exists(wDir))
    	{
            print(paste0(wDir," exists!"))
            next
    	} else {
            #if(!dir.exists(wDir)){dir.create(wDir, mode="7777", recursive=TRUE);}
            dir.create(wDir, mode="7777", recursive=TRUE)
        }
    	setwd(wDir)    

	peakFile <- paste0(mainDir,"/loopFiles/",peakFile)

    	sourcePeaks1 <- read.table(peakFile)[1:3]
    	colnames(sourcePeaks1) <- c('chrom','start','end')
    	print(head(sourcePeaks1))
    	### Load 2nd anchors ###
    	sourcePeaks2 <- read.table(peakFile)[4:6]
    	colnames(sourcePeaks2) <- c('chrom','start','end')
    	print(head(sourcePeaks2))
    	### Load the anchors ###

    	datasets1 <- list(peaks=center1D(sourcePeaks1)[,1:3])
    	datasets2 <- list(peaks=center1D(sourcePeaks2)[,1:3])
    	head(datasets1)
    	head(datasets2)
    	cluster  <- c(TRUE,FALSE)[1]

    	# In this case we have just one peakset 'peaks'
    	for(peakSet in names(datasets1))
    	{
	    print(peakSet)
            # Generate the heatmaps for 1st anchors
	    peaks1 <- datasets1[[peakSet]]
	    col <- cols[[peakSet]]
	    print(col)
	    #quit()
	    #print(nrow(peaks1))
	    #print(ncol(peaks1))
	    singlePeakDataFile1 <- 'spData1.bed'
	    # Initialize the table with the peaks1
	    write.table(center1D(peaks1[,1:3],1)[,1:3], file=singlePeakDataFile1, quote=F, row.names=F, col.names=F, sep="\t")
	    # Load the peaks1. PlotSetArray is a function from the SeqPlots package which allows to load the bigWig files with the Chip-seq tracks (first argument: tracks) 
	    # and associate them to the peakset (second argument: features), which are the peaks in this case, on the reference (refgenome), which is 'dm6' in this case.
	    # It returns the heatmaps in sqData1
	    sqData1 <- getPlotSetArray(paste0(bwDir,testBW), singlePeakDataFile1, refgenome=assembly, bin = binSize, rm0 = TRUE, ignore_strand = TRUE, xmin = range, xmax = range, type = "mf", add_heatmap = TRUE)

	    ### Generate the heatmaps for 2nd anchors ###
	    peaks2 <- datasets2[[peakSet]]
	    #col <- cols[[peakSet]]
	    #print(nrow(peaks2))
	    #print(ncol(peaks2))
	    singlePeakDataFile2 <- 'spData2.bed'
	    # Initialize the table with the peaks2
	    write.table(center1D(peaks2[,1:3],1)[,1:3], file=singlePeakDataFile2, quote=F, row.names=F, col.names=F, sep="\t")
	    # Load the peaks2. PlotSetArray is a function from the SeqPlots package which allows to load the bigWig files with the Chip-seq tracks (first argument: tracks) 
	    # and associate them to the peakset (second argument: features), which are the peaks in this case, on the reference (refgenome), which is 'dm6' in this case.
	    # It returns the heatmaps in sqData2
	    sqData2 <- getPlotSetArray(paste0(bwDir,testBW), singlePeakDataFile2, refgenome=assembly, bin = binSize, rm0 = TRUE, ignore_strand = TRUE, xmin = range, xmax = range, type = "mf", add_heatmap = TRUE)
	    ### Generate the heatmaps for 2nd anchors ###

	    ### Preparing the heatmaps for plotting for 1st and 2nd anchors data for clustering based on 1st anchors ###
	    print("Read the BW one-by one and prepare the data in allMeanData for the clustering based on anchors")
	    # allMeanData will have the targets as columns and the peaks as rows
	    allData1 <- list() # List to store the heatmaps of 1st anchors
	    allData2 <- list() # List to store the heatmaps of 2nd anchors
	    allMeanData <- c() # Vector to store the mean controbution for clustering based on 1st anchors (no 1 or 2 index needed)
	    for (set in testBW)
	    {
		print(paste0("Analysing bigWig ",set,sep=""))

	    	# Set the name of the column to plot
	    	target   <- unlist(strsplit(set,'_'))[1:3]
	    	sName <- paste(target,collapse='_')
	    	sName <- gsub('.bigWig','',sName)
	    	sName <- gsub('_merge_dm6_GS','',sName)
	    	sName <- gsub('_Embryo','',sName)	    	    
	    	sName <- gsub('.bw','',sName)
	    	print(paste0("Plot column ID: ",sName,sep=""))	    
	    	#quit()

	    	# Get the heatmap for 1st anchors
	    	# The following is a table for the bigWig set. It contains:
	    	# $anno    -> Genomic location of the peak : chr centre
	    	# $heatmap -> heatmap per peak of the target bigWig set
	    	hData1 <- sqData1$data$spData[[which(testBW == set)]]$heatmap # hData1 is the heatmap for target set. It has Nrows=Npeaks and Ncolumns=Nbins
	    	hData1[is.na(hData1)] <- 0          # Substitute nans with zeros
	    	hData1 <- clipMatrix(hData1,0.999) # Don't consider outliers

	    	# Get the heatmap for 2nd anchors
	    	hData2 <- sqData2$data$spData[[which(testBW == set)]]$heatmap # hData2 is the heatmap for target set. It has Nrows=Npeaks and Ncolumns=Nbins
	    	hData2[is.na(hData2)] <- 0          # Substitute nans with zeros
	    	hData2 <- clipMatrix(hData2,0.999) # Don't consider outliers

	    	# Select data for clustering analysis (only for 1st anchors)
	    	# cData (clusteringData) contains per each peak 1->Npeaks the elements the contribution of the target (in the bigWig)
	    	# around the peak centre (ncol(hData)/2)) from -meanBins to +minBins
	    	# summed (rowSums) and log2-transformed. So, allMeanData will contain for each peak one value per target! 
	    	if(ClusteringAnchor == 1)
	    	{
		    cData <- log2(rowSums(hData1[,seq(round(ncol(hData1)/2)-meanBins,round(ncol(hData1)/2)+meanBins)])+1)
	    	} 
	    	if(ClusteringAnchor == 2)
	    	{
                    cData <- log2(rowSums(hData2[,seq(round(ncol(hData2)/2)-meanBins,round(ncol(hData2)/2)+meanBins)])+1)
	    	}
	    	if(ClusteringAnchor == 0)
	    	{
                    cData <- cbind(log2(rowSums(hData1[,seq(round(ncol(hData1)/2)-meanBins,round(ncol(hData1)/2)+meanBins)])+1),log2(rowSums(hData2[,seq(round(ncol(hData2)/2)-meanBins,round(ncol(hData2)/2)+meanBins)])+1))
	    	    colnames(cData) <- c(paste0(sName,"_1"),paste0(sName,"_2"))
	        }	
	    	if(is.null(allMeanData))
	    	{
                    #if allMeanData is empty initilise it
		    if(ClusteringAnchor == 0)
		    {
	                #allMeanData <- matrix(cData, ncol=2)
		    	allMeanData <- cData
		    } else {
                        allMeanData <- matrix(cData, ncol=1)
		    }
	    	}else{
		    #else add a new column (cbind) to allMeanData
		    if(ClusteringAnchor == 0)
                    {
	                allMeanData <- cbind(allMeanData,cData)
		    } else {
       	                allMeanData <- cbind(allMeanData,matrix(cData, ncol=1))
		    }
	        }
		if(ClusteringAnchor == 0)
	    	{
	            print(ncol(allMeanData))
		    print(ncol(allMeanData)-1)
		    print(colnames(allMeanData))
	    	} else {
	            colnames(allMeanData)[ncol(allMeanData)] <- sName
	    	}

	    	# Store in allData1 and allData2 the heatmaps for 1st and 2nd anchors
	    	print(dim(hData1))
	    	allData1[[sName]] <- hData1
	    	print(dim(hData2))
	    	allData2[[sName]] <- hData2
  	    }
	    print(paste0("Dimension of allData of 1st anchors ",dim(allData1)))
	    print(paste0("Dimension of allData of 2nd anchors ",dim(allData2)))
	    #print(head(cData))
	    #print(tail(cData))
	    allMeanData[is.na(allMeanData)] <- 0 # Substitute nans with zeros in allMeanData.
	print(head(allMeanData))	
	### Preparing the heatmaps for plotting for 1st and 2nd anchors data for clustering based on 1st anchors ###

	# layout divides the device (R 'terminal' for a figure) into as many rows and columns 
	# as there are in matrix mat, with the column-widths and the row-heights specified in the respective arguments.
	layout(matrix(1:length(testBW),ncol=length(testBW)))

   	### Performing clustering based on 1st anchors ###
	if(peakSet == names(datasets1)[1])
	{
            # Since dist() computes the distance between rows, so transposing allMeanData, we cluster the targets (rows in the transposed) depending on the contribution per each peak (columns in the transposed) using a hierarchical clustering
	    clOrder <- hclust(dist(t(allMeanData)))$order
	    print(paste0("Cluster order ",clOrder))
		
        # Using allMeanData as it is we now cluster the peaks (rows) depending on the contribution per target (columns) using a kmeans clustering.
	if(cluster == TRUE)
        {
	    set.seed(32)
            # kmeans cluster the peaks in kMax clusters using a kmeans algorithm depending on the value of each peak for the different targets
	    #	    kMeans <- kmeans(allMeanData, kMax)
	    # Test the best number of clusters
	    distances = dist(allMeanData)
	    kMeans   <- hclust(distances, method="ward.D2")
	    calinhara_scores <- foreach(n=1:kMax) %dopar% {	    
	    	print(paste0("Computing CH score of the clustering in ",n," clusters"))
		#kMeans <- kmeans(allMeanData, n)
		clustering <- cutree(kMeans, k = n)
		#clustering <- kMeans$cluster
		score      <- calinhara(distances, clustering, cn=n)
		if(is.na(score)){score=0;}
		print(paste0("CH score of the clustering in ",n," clusters ",score))	    
		score
	    }
	    print(calinhara_scores)
	    score <- unlist(calinhara_scores) # Variable with the the Calinhara scores per each number of clusters
	    print(max(score))
	    print(score)	
	    #best_nclusters = which(score == max(score))
	    best_nclusters = 2 #which(score == max(score))	    
	    cat(paste("#Best number of clusters ",best_nclusters), "\n")
	    kMeans <- kmeans(allMeanData, best_nclusters)
	    # kMeans$cluster <- cutree(kMeans, k = best_nclusters)
	    # kOrder is the order of the peaks after the kmeans clustering
            # kMeans$cluster is a vector of integers (from 1:kMax) indicating the cluster to which each point is allocated.
            # order() orders the peaks of allMeanData depending on the cluster they belong.
            # That is: first all the peaks of cluster 1, then the ones of cluster 2 and so on...
            kOrder <- order(kMeans$cluster,apply(allMeanData,1,max))
            print(paste0("Order ",kOrder))
        }else{
	    print('NOT CLUSTERING')
	    kMeans <- km[[peakSet]]
	    kOrder <- order(kMeans$cluster,apply(allMeanData,1,max))
	}
        peaks1$k           <- kMeans$cluster
	peaks2$k           <- kMeans$cluster
        
        # Number of heatmaps-panels (columns) in the final figure = (number of BigWigs) * (number of peak-sets) * (2 loop anchors)
        nPanels <- length(testBW)*length(names(datasets1))*2 + 1 + 1
        # Width of the output figure
        outWidth <- nPanels*1.32
 	
        # Open the R device on which to save the figure as a .pdf
	#pdf(file=paste(paste('SeqPlots',paste('meanBins',meanBins,sep=""),paste('Nclusters',best_nclusters,sep=""),sep='_'),'pdf',sep="."),width=outWidth,height=9)
	pdf(file=paste0('SeqPlots_meanBins',meanBins,'_Nclusters',best_nclusters,'.pdf'),width=outWidth,height=9)	
		
	layout(matrix(1:(nPanels*2),ncol=nPanels), heights=c(outWidth/3.75,1))
        #par(mai=rep(0.125,4))
	#par sets the margins (mai=c(bottom, left, top, right)) and the stiles of y and x axis (Style "i" (internal) just finds an axis with pretty labels that fits within the original data range)
	par(mai=rep(c(0.3,0.125),2), yaxs="i", xaxs="i")
    }
    ### Performing clustering based on 1st anchors ###

    ### Plotting the obtained cluster analysis including also 2nd anchors ###
    # Plot the leftmost column with the labels of the clusters
    nameCex <- 1.25
    plot(0,xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), col='white', frame=F)
    clLabels <- sapply(max(kMeans$cluster):1,function(k){
	      y <- (max(which(sort(kMeans$cluster, decreasing=FALSE) == k)) - length(which(kMeans$cluster == k))/2)/length(kMeans$cluster)
       	      txt <- paste0(paste0('CL',k),'\n',length(which(kMeans$cluster == k)),' loops')
              text(0.5,y,txt,cex=nameCex)
	      return(txt)})
    addClusterLines(kMeans,FALSE,FALSE,col=1,max=length(kMeans$cluster))
    plot(0,xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), col='white', frame=F)
    
    print("Plot the columns of the figure one per target in the order given by hclust (clOrder)")
    print(clOrder)
    if(ClusteringAnchor == 0)
    {
        clOrder = clOrder[clOrder<=length(testBW)]
    }
    print(clOrder)    
    
    for(set in names(allData1)[clOrder])
    {
        #title <- gsub(pattern=c(assembly,'H3'),replacement=c('',''),set)	
	#text(0.5,0.0,title,cex=nameCex)

    	### Plots for 1st anchors ###
	# Add heatmaps of the clustered peaks for 1st anchors
	# matrix1 is the heatmap for the target (e.g., histone mark H4K4me3) named set for 1st anchors
        matrix1 <- allData1[[set]]
        #print(dim(matrix1))
        # Use clipMatrix to remove outliers    
        matrix1 <- clipMatrix(matrix1,0.99)
        # Order the rows of matrix1 (one per peak) depending on the kmeans-cluster they belong and take the transpose matrix
        matrix1 <- t(matrix1[kOrder,])
        print(paste("Dimension of the heatmap ",dim(matrix1)[[1]],dim(matrix1)[[2]],sep=" "))       
        # Set the title of the column
	#title <- "1st anchor"
	title <- paste0(gsub(pattern=c("_dm6"),replacement='',set),"_1")
        # Plot matrix matrix1 using color-scale col
	image(matrix1, col=col, useRaster=raster, axes=F, main=title, cex.main=1.5)

        # Add black lines to mark the cluster borders
	addClusterLines(kMeans,FALSE,FALSE,col=1,max=length(kMeans$cluster))
        # Add yellow lines to mark the bins around the centre of the peak used for clustering
 	if(ClusteringAnchor == 1)
        {
            abline(v=(nrow(matrix1)/2-meanBins)/nrow(matrix1),lty=2,col='yellow',lwd=0.5)
 	    abline(v=(nrow(matrix1)/2+meanBins)/nrow(matrix1),lty=2,col='yellow',lwd=0.5)
        }
 	if(ClusteringAnchor == 0)
        {
            abline(v=(nrow(matrix1)/2-meanBins)/nrow(matrix1),lty=2,col='yellow',lwd=0.5)
 	    abline(v=(nrow(matrix1)/2+meanBins)/nrow(matrix1),lty=2,col='yellow',lwd=0.5)
        }	
	if(cluster == TRUE)
        {
	    clOrd <- sort(kMeans$cluster)
	}else{
	    kMeans <- km[[peakSet]]
	    clOrd <- kMeans$cluster[kOrder]
	}

        print(title)
	for(x in seq(from=1,to=max(kMeans$cluster),by=1))
	{
 	    
	}
	print(nrow(matrix1[,which(clOrd == 1)]))
	print(nrow(matrix1[,which(clOrd == 2)]))
	print(nrow(matrix1[,which(clOrd == 3)]))
	print(nrow(matrix1[,which(clOrd == 4)]))
	print(nrow(matrix1[,which(clOrd == 5)]))	
        #yMax <- max(unlist(sapply(1:max(kMeans$cluster),function(x){max(smooth(rowMeans(matrix1[,which(clOrd == x)],na.rm=T)))})))
        #yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(rowMeans(matrix1[,which(clOrd == x)],na.rm=T)))})))
        yMax <- max(unlist(sapply(1:max(kMeans$cluster),function(x){max(smooth(safeRowMeans(matrix1[,which(clOrd == x), drop=FALSE])))})))
	yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(safeRowMeans(matrix1[,which(clOrd == x), drop=FALSE])))})))

        for (i in 1:max(kMeans$cluster))
        {
	    print(paste0("Cluster ",i))
	    clMtx <- matrix1[,which(clOrd == i), drop=FALSE]
	    if(i==1)
            {
	    	#plot(smooth(rowMeans(clMtx,na.rm=T)), type='l', lwd=2, xlab='', ylab='', xaxt='n', yaxt='n',col=i,ylim=c(yMin,yMax),main=gsub('H3','',set))
		plot(smooth(safeRowMeans(clMtx)), type='l', lwd=2, xlab='', ylab='', xaxt='n', yaxt='n',col=i,ylim=c(yMin,yMax),main=gsub('H3','',set))		
		axis(1,labels=c(paste0('-',range/1e3,'kb'),0,paste0(range/1e3,'kb')), at=c(1,nrow(clMtx)/2,nrow(clMtx)))
	    }else{
                #lines(smooth(rowMeans(clMtx,na.rm=T)), lwd=2,col=i)
    		lines(smooth(safeRowMeans(clMtx)), lwd=2,col=i)		
	    }
	}
	### Plots for 1st anchors ###

	### Plots for 2nd anchors ###
        # Add heatmaps of the clustered peaks for 2nd anchors
        # matrix2 is the heatmap for the target (e.g., histone mark H4K4me3) named set for 2nd anchors
        matrix2 <- allData2[[set]]
        print(dim(matrix2))
    	# Use clipMatrix to remove outliers    
	matrix2 <- clipMatrix(matrix2,0.99)
	# Order the rows of matrix2 (one per peak) depending on the kmeans-cluster they belong and take the transpose matrix
    	matrix2 <- t(matrix2[kOrder,])
    	print(paste("Dimension of the heatmap ",dim(matrix2)[[1]],dim(matrix2)[[2]],sep=" "))       
    	# Set the title of the column
	#	title <- "2nd anchor"
	title <- paste0(gsub(pattern=c('H3'),replacement='',set),"_2")
        # Plot matrix matrix2 using color-scale col
    	image(matrix2, col=col, useRaster=raster, axes=F, main=title, cex.main=1.5)

	# Add black lines to mark the cluster borders
	addClusterLines(kMeans,FALSE,FALSE,col=1,max=length(kMeans$cluster))
        # Add yellow lines to mark the bins around the centre of the peak used for clustering
	if(ClusteringAnchor == 2)
        {
	    abline(v=(nrow(matrix2)/2-meanBins)/nrow(matrix2),lty=2,col='yellow',lwd=0.5)
 	    abline(v=(nrow(matrix2)/2+meanBins)/nrow(matrix2),lty=2,col='yellow',lwd=0.5)
        }
	if(ClusteringAnchor == 0)
        {
	    abline(v=(nrow(matrix2)/2-meanBins)/nrow(matrix2),lty=2,col='yellow',lwd=0.5)
 	    abline(v=(nrow(matrix2)/2+meanBins)/nrow(matrix2),lty=2,col='yellow',lwd=0.5)
        }	
	if(cluster == TRUE)
        {
	    clOrd <- sort(kMeans$cluster)
	}else{
	    kMeans <- km[[peakSet]]
	    clOrd <- kMeans$cluster[kOrder]
	}

        print(title)
        #yMax <- max(unlist(sapply(1:max(kMeans$cluster),function(x){max(smooth(rowMeans(matrix2[,which(clOrd == x)],na.rm=T)))})))
        #yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(rowMeans(matrix2[,which(clOrd == x)],na.rm=T)))})))
        yMax <- max(unlist(sapply(1:max(kMeans$cluster),function(x){max(smooth(safeRowMeans(matrix2[,which(clOrd == x), drop=FALSE])))})))
        yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(safeRowMeans(matrix2[,which(clOrd == x), drop=FALSE])))})))	

        for (i in 1:max(kMeans$cluster))
        {
	    print(paste0("Cluster ",i))
	    clMtx <- matrix2[,which(clOrd == i), drop=FALSE]
	    if(i==1)
            {
		#plot(smooth(rowMeans(clMtx,na.rm=T)), type='l', lwd=2, xlab='', ylab='', xaxt='n', yaxt='n',col=i,ylim=c(yMin,yMax),main=gsub('H3','',set))
		plot(smooth(safeRowMeans(clMtx)), type='l', lwd=2, xlab='', ylab='', xaxt='n', yaxt='n',col=i,ylim=c(yMin,yMax),main=gsub('H3','',set))		
		axis(1,labels=c(paste0('-',range/1e3,'kb'),0,paste0(range/1e3,'kb')), at=c(1,nrow(clMtx)/2,nrow(clMtx)))
	    }else{
                #lines(smooth(rowMeans(clMtx,na.rm=T)), lwd=2,col=i)
                lines(smooth(safeRowMeans(clMtx)), lwd=2,col=i)		
	    }
	}
	### Plots for 2nd anchors ###

    }
    dev.off()

    # Write the loops per cluster
    #peaks1df = do.call(cbind, lapply(peaks1, as.data.frame))    
    #colnames(peaks1df) = c('chrom1','start1','end1','k1')
    #peaks2df = do.call(cbind, lapply(peaks2, as.data.frame))
    #colnames(peaks2df) = c('chrom2','start2','end2','k')
    loops <- cbind(peaks1,peaks2)
    colnames(loops) <- c('chrom1','start1','end1','k1','chrom2','start2','end2','k')
    print(head(loops))
	for (i in 1:(max(kMeans$cluster)))
        {
	    filename <- paste(paste("loops",condition,"in_cluster",i,"MeanBins",meanBins,"Nclusters",best_nclusters,sep="_"),"bed",sep=".")
            print(filename)
            write.table(loops[loops$k == i,c("chrom1", "start1", "end1","chrom2","start2","end2","k")],file=filename, quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE)
    	}
        setwd(paste0(mainDir,"/",condition,"/"))
    }

}
}
