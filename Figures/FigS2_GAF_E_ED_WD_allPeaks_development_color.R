# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("seqplots")
# 
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github('przemol/seqplots', build_vignettes=FALSE)
# 
# library(seqplots)

setwd("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/SeqPlots/GAF_ED_WD_E_development/")

require(seqplots)
library(seqplots)

source('/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2024_Secondary_analysis/seqPlots/auxFunctions.R')
source("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2024_Secondary_analysis/seqPlots/auxSeqplotfunction.R")


# target <- 'E_GAF_allPeaks'
target <- 'GAF_E_ED_WD_allPeaks_development'
runName <- target
data <- read.delim('/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/SeqPlots/GAF_ED_WD_E_development/E_ED_WD_merged_GAF_narrowPeaks_s100_FC2.bed', header = F)
data <- subset(data, select = c(1:3))
colnames(data)[1:3] <- c('chrom','start','end')
peaks <- data 

clType <- c('INS','EPI')[2]
nSel <- 1:2
if(clType == 'INS'){
  nSel <- 1:4
}else{
  nSel <- 1:2
}

bwDir <- '/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/ChIP_CnR_new_tracks/'
allBW <- list.files(bwDir, pattern='bigWig|bw')
testBW <- allBW
testBW <- testBW[grep('ATAC|GAF|K27|K4|PH',testBW)]
print(sort(testBW))

range <- 5e3 # 2e3
binSize <- range/100
raster <- TRUE

cName <- paste0(runName,'')

testFile <- 'spData.bed'
options(scipen=15)
write.table(peaks, file=testFile, quote=F, row.names=F, col.names=F, sep="\t")

sqData <- getPlotSetArray2(paste0(bwDir,testBW), testFile, refgenome='dm6', bin = binSize, rm0 = TRUE, ignore_strand = TRUE, xmin = range, xmax = range, type = "mf", add_heatmap = TRUE)

allData <- list()
allMeanData <- c()
fullMatrix <- c()
for (set in testBW){
  sName <- gsub('.bigWig','',paste(unlist(strsplit(set,'_'))[nSel],collapse='_'))
  hData <- sqData$data$spData[[which(testBW == set)]]$heatmap
  hData <- clipMatrix(hData,0.999)
  allData[[sName]] <- hData
  
  meanBins <- (range/binSize)/2
  print(set)
  cData <- hData[,seq(round(ncol(hData)/2)-meanBins,round(ncol(hData)/2)+meanBins)]
  cData <- scaleMatrix(cData,0,1,1e-6)
  fullMatrix <- cbind(fullMatrix,scaleMatrix(hData,0,1,1e-6))
  colnames(fullMatrix)[((ncol(fullMatrix)-ncol(hData))+1):ncol(fullMatrix)] <- paste(set,1:ncol(hData),sep='_')
  
  cData <- log2(rowSums(hData[,seq(round(ncol(hData)/2)-meanBins,round(ncol(hData)/2)+meanBins)])+1)
  allMeanData <- cbind(allMeanData,matrix(cData, ncol=1))
  colnames(allMeanData)[ncol(allMeanData)] <- sName
}

allMeanData[is.na(allMeanData)] <- 0
fullMatrix[is.na(fullMatrix)] <- 0
rownames(fullMatrix) <- 1:nrow(fullMatrix)

####################
set.seed(32)
kMaxPrimary <- 6

#################### Select datasets
selSet <- c('GAF|H3K27|ATAC','E_K27me3','PH','mcCov','ATAC','PC','K4','K27ac')[1]
kMeansDiff <- kmeans(allMeanData[,grep(selSet,colnames(allMeanData))], kMaxPrimary)
# kMeansDiff <- kmeans(fullMatrix[,grep(selSet,colnames(fullMatrix))], kMaxPrimary)
####################
kOrderDiff <- order(kMeansDiff$cluster)
peaks$k <- kMeansDiff$cluster
####################

#############
data <- peaks
#############
clOrder <- c(1:length(testBW))

if(clType == 'INS'){
  clOrder <- unlist(sapply(c('mC_Embryo','mC_LD','mC_WD'),function(x){which(startsWith(testBW,x))}))
}else{
   # clOrder <- c(unlist(sapply(c('GAF','PH','K27ac','K27me3','ATAC'),function(x){grep(x,testBW)})))
   clOrder <- c(unlist(sapply(c('E_GAF','ED_GAF','WD_GAF','E_PH','ED_PH','WD_PH','E_H3K27ac','ED_H3K27ac','WD_H3K27ac','E_H3K27me3','ED_H3K27me3','WD_H3K27me3','E_ATAC','ED_ATAC','WD_ATAC'),function(x){grep(x,testBW)})))
}

#############
#color palette
# col <- colorRampPalette(c('grey90','darkblue'))(200)
# col <- colorRampPalette(c('grey90','grey15'))(200)

track_colors <- list(
  "E_GAF"     = colorRampPalette(c("grey90", "goldenrod"))(200),
  "WD_GAF"     = colorRampPalette(c("grey90", "goldenrod"))(200),
  "ED_GAF"     = colorRampPalette(c("grey90", "goldenrod"))(200),
  "E_H3K27me3"= colorRampPalette(c("grey90", "steelblue3"))(200),
  "WD_H3K27me3"= colorRampPalette(c("grey90", "steelblue3"))(200),
  "ED_H3K27me3"= colorRampPalette(c("grey90", "steelblue3"))(200),
  "E_ATAC"    = colorRampPalette(c("grey90", "violetred"))(200),
  "WD_ATAC"    = colorRampPalette(c("grey90", "violetred"))(200),
  "ED_ATAC"    = colorRampPalette(c("grey90", "violetred"))(200),
  "E_H3K27ac"   = colorRampPalette(c("grey90", "tomato"))(200),
  "WD_H3K27ac"   = colorRampPalette(c("grey90", "tomato"))(200),
  "ED_H3K27ac"   = colorRampPalette(c("grey90", "tomato"))(200),
  "E_PH"      = colorRampPalette(c("grey90", "springgreen4"))(200),
  "ED_PH"      = colorRampPalette(c("grey90", "springgreen4"))(200),
  "WD_PH"      = colorRampPalette(c("grey90", "springgreen4"))(200)
  # Add more mappings as needed
)
####

extraColumns <- 1
mainCex <- 1.4
ph <- 13 # 12
pw <- (extraColumns+length(clOrder))*(ph/10) # (extraColumns+length(clOrder))*(ph/10)

cName <- paste0(cName,'_k_',kMaxPrimary)

pdf(file=paste0(cName,'.pdf'),width=pw,height=ph)
layout(matrix(1:((length(clOrder)+extraColumns)),ncol=length(clOrder)+extraColumns))
par(mai=rep(c(0.3,0.1),2), xaxs="i", yaxs="i", family='sans')

nameCex <- 1.5
plot(0,xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), col='white', frame=F)
clLabels <- sapply(max(kMeansDiff$cluster):1,function(k){
  y <- (max(which(sort(kMeansDiff$cluster, decreasing=FALSE) == k)) - length(which(kMeansDiff$cluster == k))/2)/length(kMeansDiff$cluster)
  
  if(length(which(kMeansDiff$cluster == k))/length(kMeansDiff$cluster) < 0.05){
    txt <- paste0(paste0('CL',k),' - ',length(which(kMeansDiff$cluster == k)))
  }else{
    txt <- paste0(paste0('CL',k),'\n',length(which(kMeansDiff$cluster == k)),'\nPEAKS')
  }
  
  text(0.5,y,txt,cex=nameCex)
  return(txt)
})
addClusterLines(kMeansDiff,FALSE,FALSE,col=1,max=length(kMeansDiff$cluster))

for (set in names(allData)[clOrder]){
  
  m <- allData[[set]]
  m[is.na(m)] <- 0
  m <- clipMatrix(m,0.99)
  m <- t(m[kOrderDiff,])
  
  title <- unlist(strsplit(set,'_ENC'))[1]
  title <- gsub('_NA|_INS|mC_','',title)
  # image(m, col=col, useRaster=raster, axes=F, main=gsub('H2A|H3','',title), cex.main=mainCex) #, zlim=zLim)
  # Extract a key from the set name to choose the color
  track_key <- names(track_colors)[sapply(names(track_colors), function(x) grepl(x, set))][1]
  
  # Fallback if no match is found
  track_palette <- if (!is.na(track_key)) track_colors[[track_key]] else colorRampPalette(c("grey95", "black"))(200)
  
  image(m, col=track_palette, useRaster=raster, axes=F, main=gsub('','',title), cex.main=mainCex)
  
  addClusterLines(kMeansDiff,FALSE)
  if(clType == 'INS'){
    if(length(grep('INS',set)) > 0){
      abline(v=(nrow(m)/2-meanBins)/nrow(m),lty=2,col='yellow',lwd=0.5)
      abline(v=(nrow(m)/2+meanBins)/nrow(m),lty=2,col='yellow',lwd=0.5)
      axis(1,paste0((2*meanBins*binSize)/1e3,'kb'),at=0.5,tick=FALSE)
    }
  }else{
    if(set == names(allData)[clOrder][1]){
      axis(1,c(paste0(-range/1e3,'kb'),paste0(range/1e3,'kb')),at=c(0,1),las=1)
    }
  }
}
dev.off()


next

# Write the peaks per cluster
write.table(peaks,file=paste0(cName,'.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==1,], file=paste0(cName,'_cl1.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==2,], file=paste0(cName,'_cl2.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==3,], file=paste0(cName,'_cl3.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==4,], file=paste0(cName,'_cl4.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==5,], file=paste0(cName,'_cl5.bed'), quote=F, row.names=F, col.names=F, sep="\t")
write.table(peaks[peaks$k ==6,], file=paste0(cName,'_cl6.bed'), quote=F, row.names=F, col.names=F, sep="\t")


# ########## get the average profiles
# ext <- range
# clip <- 0.995
# 
# avgSets <- c("E_GAF","ED_GAF","WD_GAF","E_PH","ED_PH","WD_PH")
# 
# # Clip intensities but don't rescale yet
# for(set in avgSets){
#   allData[[set]] <- clipMatrix(allData[[set]], clip)
# }
# 
# # Compute a global max to scale all together
# global_max <- max(sapply(allData[avgSets], max, na.rm = TRUE))
# 
# for(set in avgSets){
#   allData[[set]] <- allData[[set]] / global_max
# }
# 
# 
# pdf(file=paste0(cName,'_avgProf_GAF_PH.pdf'),width=ph/4,height=ph)
# layout(matrix(1:(max(kMeansDiff$cluster)+1),ncol=1))
# 
# par(mai=rep(0.3,4))
# 
# ##change reference labels
# plot(0,col='white',xlab='', ylab='', yaxt='n', xaxt='n', frame=F, ylim=c(0,5), xlim=c(0,5))
# text(avgSets[1],col='yellow2',x=2.5,y=0, pos=3, cex=1.6)
# text(avgSets[2],col='goldenrod2',x=2.5,y=1, pos=3, cex=1.6)
# text(avgSets[3],col='goldenrod4',x=2.5,y=1, pos=3, cex=1.6)
# text(avgSets[4],col='chartreuse',x=2.5,y=2, pos=3, cex=1.6)
# text(avgSets[5],col='darkolivegreen3',x=2.5,y=2, pos=3, cex=1.6)
# text(avgSets[6],col='springgreen4',x=2.5,y=3, pos=3, cex=1.6)
# 
# 
# for(k in max(kMeansDiff$cluster):1){
# 
#   for (set in avgSets){
#     selRows <- which(kMeansDiff$cluster == k)
#     m <- allData[[set]]
#     m <- m[selRows,]
#     m <- t(m)
# 
#     lwd <- 2.5
# 
#     # y <- rowMeans(m,na.rm=T)
#     # ylim <- c(0,1)
#     y <- rowMeans(m, na.rm = TRUE)
#     ylim <- c(0, max(y) * 2.4)  # 
# 
#     col <- 'yellow2'
#     if(set == avgSets[2]){col <- 'goldenrod2'}
#     if(set == avgSets[3]){col <- 'goldenrod4'}
#     if(set == avgSets[4]){col <- 'chartreuse'}
#     if(set == avgSets[5]){col <- 'darkolivegreen3'}
#     if(set == avgSets[6]){col <- 'springgreen4'}
#     if(set == avgSets[1]){
#       plot(y,type='l', lwd=1.5, col=col,cex.main=1.4,xlab='', ylab='', yaxt='n', xaxt='n', frame=F, ylim=ylim, main=paste0('CL ',k))
#       axis(1,labels=c(-ext,0,ext),at=c(0,nrow(m)/2,nrow(m)),cex=1.2)
#       abline(v=nrow(m)/2,col='darkslategrey',lty=2,lwd=2)
#       lines(y,col=col,lwd=lwd)
#     }else{
#       lines(y,col=col,lwd=lwd)
#     }
#   }
# }
# 
# dev.off()
