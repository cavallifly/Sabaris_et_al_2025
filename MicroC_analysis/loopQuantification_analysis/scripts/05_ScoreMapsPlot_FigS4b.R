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
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

##############################################################################
geneCoordinates <- gintervals.load("intervals.ucscCanGenes")
rownames(geneCoordinates) <- geneCoordinates$geneName
##############################################################################

k    <- 250
kexp <- 2*k
zmin <- -100
zmax <-  100
scoreData <- list(
        larvae_DWT       = paste0("hic.hicScores_larvae_DWT_merge_dm6_BeS.hicScores_larvae_DWT_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb"),
	LD_DPRE2         = paste0("hic.hicScores_LD_DPRE2_merge_dm6_BeS.hicScores_LD_DPRE2_merge_dm6_BeS_score_k",k,"_kexp",kexp,"_step1Mb"),
        )

extentions  <- c(1,2)

for(resolutions in list(c(3.0e3,3.0e3)))
{
    print(resolutions)

    #for(genomicLocus in c("dac","en","NetA","Bx-C"))
    for(genomicLocus in c("dac"))    
    {

        if(genomicLocus == "dac")
	{
	    # Genomic locations of dac locus
	    chr <- "chr2L"
	    start1 <- 16000000
	    end1   <- 17000000
	    start2 <- 16321000
	    end2   <- 16512000
	    pre1   = gintervals("chr2L",16422435,16423525)
	    pre2   = gintervals("chr2L",16485929,16486572)
	    boundary1 <- gintervals("chr2L",16354000,16354001)
	    boundary2 <- gintervals("chr2L",16494000,16494001)
	    geneAnnotation    <- g2d(rbind(geneCoordinates["dac",1:3],expand1D(gintervals("chr2L",16422514,16422515),1e3),expand1D(gintervals("chr2L",16560374,16560375),1e3))[,1:3])
	}

	if(genomicLocus == "en")
        {
	    # Genomic locations of en locus
	    chr <- "chr2R"
	    start1 <- 11000000
	    end1   <- 12000000
	    start2 <- 11421500
	    end2   <- 11612500
	    pre1   = gintervals("chr2R",11476251,11476252)
	    pre2   = gintervals("chr2R",11528751,11528752)
	    boundary1 <- gintervals("chr2R",11452000,11452001)
	    boundary2 <- gintervals("chr2R",11582000,11582001)
	    #2R:11,524,006..11,528,210 [-]
	    geneAnnotation   <- g2d(rbind(geneCoordinates["en",1:3],expand1D(gintervals("chr2R",11524006,11524007),1e3),expand1D(gintervals("chr2R",11528210,11528211),1e3))[,1:3])
	}

	if(genomicLocus == "NetA")
        {
	    # Genomic locations of NetA locus
	    ## NetA locus: (you already did for the reviewer figure of NSMB and GB) 
	    ## ChrX: 14 600 000..14 800 000
	    chr <- "chrX"
	    start1 <- 14600000
	    end1   <- 14800000
	    start2 <- 14600000
	    end2   <- 14800000
	    pre1   = gintervals("chrX",14656251,14656252)
	    pre2   = gintervals("chrX",14751251,14751252)
	}

	if(genomicLocus == "Bx-C")
        {
	    # Genomic locations of Bx-C locus
	    ## Bx-C:
	    ## 3R: 16 625 000..17 000 000
	    ## Indicate Position of Fab7 PRE: 16 899 628..16 900 158
	    chr <- "chr3R"
	    start1 <- 16625000
	    end1   <- 17000000
	    start2 <- 16625000
	    end2   <- 17000000
	    pre1   = gintervals("chr3R",16899628,16899629)
            pre2   = gintervals("chr3R",16900158,16900159)
	}

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

		outFile = paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_FigS4b_annot_",annotated,".png")
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
	    	    outMatrix <- paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_FigS4b_",set,".tsv")
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

		png(file=paste0("scoreMapsk",k,"kexp",kexp,"_",chr,"_",start,"bp_",end,"bp_r",res,"bp_FigS4b_annot_",annotated,"_toModify.png"),width=nc*slotSize,height=nr*slotSize)

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
pdf(paste0('scores_colorbar_FigS4b.pdf'),width=8,height=8)
labels <- c(as.character(zmin),as.character(as.integer(zmin/2)),"0",as.character(as.integer(zmax/2)),as.character(zmax))
x <- length(scoreCol)
plot(y=1:x,x=rep(0,x),col=scoreCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
dev.off()
warnings()
