library(misha)
library(shaman)
library(doParallel)

##########################################
### misha working DB
mDBloc <-  "/zdata/data/mishaDB/trackdb/"
db <- "dm6"
dbDir <- paste0(mDBloc,db,"/")
gdb.init(dbDir)
gdb.reload()

source("/zdata/data/auxFunctions/auxFunctions.R")
options(scipen=20,gmax.data.size=1.e8,shaman.sge_support=1)

##############################################################################
geneCoordinates <- gintervals.load("intervals.ucscCanGenes")
rownames(geneCoordinates) <- geneCoordinates$geneName
##############################################################################

k    <- 250
kexp <- 2*k
zmin <- -100
zmax <-  100
scoreData <- list(
	ED_PH18       = paste0("microC.mCScores_ED_PH18_merge_dm6_BeS.mCScores_ED_PH18_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb"),
	ED_PH29       = paste0("microC.mCScores_ED_PH29_merge_dm6_BeS.mCScores_ED_PH29_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb")
        )

extentions  <- c(1,2)

inRegions = "Furlong_allLoops_regions.tab"
genomicLoci = read.table(inRegions, header=F)
colnames(genomicLoci) <- c("chrom","start","end","region","genes")
print(head(genomicLoci))

inTADs    = "tad_boundaries.tsv"
TADboundaries = read.table(inTADs, header=T)
print(head(TADboundaries))

for(resolutions in list(c(1.0e3,1.0e3),c(3.0e3,3.0e3)))
{
    print(resolutions)

    for(r in 1:nrow(genomicLoci))
    {
        genomicLocus <- genomicLoci[r,]
	print(genomicLocus)

	chr    <- genomicLocus$chrom
	chrSize <- gintervals.all()[gintervals.all()$chrom == chr,]$end
	print(paste0("Chromosome ",chr," size ",chrSize))	

	region <- genomicLocus$region
	genes  <- genomicLocus$genes
	print(paste0(genomicLocus$chrom," ",genomicLocus$start," ",genomicLocus$end))

	start1 <- genomicLocus$start - 250000
	end1   <- genomicLocus$end   + 250000
	if(start1 < 0)
	{
	    start1 <- 10000
	    end1   <- genomicLocus$end + genomicLocus$start
	}
	if(end1 > chrSize)
	{
	    end1   <- as.integer(chrSize/resolutions[[1]])*resolutions[[1]] - resolutions[[1]]
	    start1 <- end1 - 500000
	}

	start2 <- start1
        end2   <- end1
	print(paste0(genomicLocus$chrom," ",genomicLocus$start," ",genomicLocus$end," ",start1," ",end1," ",start2," ",end2), quote=F)

        pre1   = gintervals(genomicLocus$chrom, genomicLocus$start, genomicLocus$start+1)
        pre2   = gintervals(genomicLocus$chrom, genomicLocus$end  , genomicLocus$end+1  )
	print(pre1)
	print(pre2)	

	boundary1 = gintervals(chr, start1, start1+1)
	boundary2 = gintervals(chr, end1  , end1+1  )
	print(boundary1)
	print(boundary2)	
	
	### valid for any locus ###
	preAnnotation    <- g2d(rbind(expand1D(pre1,1e3),expand1D(pre2,1e3))[,1:3])
	boundary1Annotation <- g2d(rbind(expand1D(boundary1,1e3)))
	boundary2Annotation <- g2d(rbind(expand1D(boundary2,1e3)))
	
	for(annotated in c("FALSE","TRUE"))
	{
    	    for(extend in extentions)
    	    {
		i <- which(extentions == extend)
		res <- resolutions[i]
		if (end1 - start1 > 5000000)
		{
	    	    res <- 10000
		}

		if(extend == 1)
		{
	    	    start <- start1
	    	    end   <- end1
		}
		if(extend == 2)
		{
	    	    start <- start2
	    	    end   <- end2
		}

		intrv <- gintervals(chr,start,end)
		exIntrv <- giterator.intervals(intervals=g2d(intrv),iterator=c(res,res))

		target <- ""
		color  <- ""
		scoreCol <- bonevScale()

		### Plot Score maps
		nr <- 1
		slotSize <- 800

		nc <- length(names(scoreData))

		outFile = paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_",region,"_",genes,"_Fig4_annot_",annotated,".png")
		if(file.exists(outFile)){next;}

		png(file=outFile,width=nc*slotSize,height=nr*slotSize)	
		print(paste0("Plotting map between ",start,"-",end))

		layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
		par(mai=rep(0.45,4))

		plotData <- list()

		for(set in names(scoreData))
		{

		    print(paste0("Analysing sample ",set))
		    track  <- scoreData[[set]]
		    data   <- gextract(track,exIntrv,iterator=exIntrv,colnames="score")
	    	    data[is.na(data)] <- 0
	    	    outMatrix <- paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_",region,"_",genes,"_Fig4_",set,".tsv")
		    #if(file.exists(outMatrix)){next;}
	    	    write.table(data[,c(1,2,3,4,5,6,7)], file = outMatrix, sep="\t", row.names=FALSE, quote=FALSE)
		    
	    	    mtx    <- acast(data,start1~start2,value.var="score")

	    	    counts <- gextract(track,g2d(intrv),colnames="score")
	    	    regReads <- nrow(counts)/2

	    	    title <- paste0(set," - ", signif(regReads/1000,3),"k in region (",signif(min(data$score),2),":",signif(max(data$score),2),")")

	    	    image(mtx, main=title, useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=scoreCol, zlim=c(zmin,zmax))
	    	    axis(1,at=c(0,0.5,1),labels=c(signif(start/1e6,3),chr,signif(end/1e6,3)),cex.axis=2.6)
	    	    if (annotated == "TRUE")
	    	    {
                        # Annotation on the diagonal
			plotRectangles2D(boundary1Annotation   ,g2d(intrv),lwd=6,image=TRUE,col='black')				
			plotRectangles2D(preAnnotation   ,g2d(intrv),lwd=6,image=TRUE,col='yellow')
			plotRectangles2D(boundary2Annotation   ,g2d(intrv),lwd=6,image=TRUE,col='black')				

			# PRE loop annotation
			PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
			plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)
		    }
	        }
		dev.off()

		png(file=paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_",region,"_",genes,"_Fig4_annot_",annotated,"_toModify.png"),width=nc*slotSize,height=nr*slotSize)

		layout(matrix(1:(nr*nc),ncol=nc,byrow=TRUE),respect=TRUE)
		par(mai=rep(0.45,4))

		plotData <- list()

		for(set in names(scoreData))
		{
	    	    print(paste0("Analysing sample ",set))
		    track  <- scoreData[[set]]
	    	    data   <- gextract(track,exIntrv,iterator=exIntrv,colnames="score")
	    	    data[is.na(data)] <- 0
	    	    mtx    <- acast(data,start1~start2,value.var="score")

	    	    title <- ""

		    image(mtx, main=title, useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=scoreCol, zlim=c(zmin,zmax))
		    if (annotated == "TRUE")
		    {
                        # Annotation on the diagonal
			plotRectangles2D(boundary1Annotation   ,g2d(intrv),lwd=6,image=TRUE,col='black')				
			plotRectangles2D(preAnnotation   ,g2d(intrv),lwd=6,image=TRUE,col='yellow')
			plotRectangles2D(boundary2Annotation   ,g2d(intrv),lwd=6,image=TRUE,col='black')						

			# PRE loop annotation
			PREloop <- expand2D(gintervals.2d(pre1$chrom,pre1[2],pre1[3],pre2$chrom,pre2[2],pre2[3]),1e3)
			plotRectangles2D(PREloop,g2d(intrv),lwd=6,image=TRUE)
 	    	    }
		}
		dev.off()
	    }
	}
    }
}
# Plot the colorbar once for all
pdf(paste0('scores_colorbar_Fig4.pdf'),width=8,height=8)
labels <- c(as.character(zmin),as.character(as.integer(zmin/2)),"0",as.character(as.integer(zmax/2)),as.character(zmax))
x <- length(scoreCol)
plot(y=1:x,x=rep(0,x),col=scoreCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
dev.off()
warnings()
