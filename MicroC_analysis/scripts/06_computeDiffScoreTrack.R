require(misha)
require(shaman)

require(doParallel)

##########################################
### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('./scripts/auxFunctions.R')
options(scipen=20,gmax.data.size=0.5e8)

refData <- list(
	hic_larvae_DWT_merge_dm6_BeS            = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1'),
	hic_larvae_double_merge_dm6_BeS         = c('hic.hic_larvae_double_Rep1_dm6_BeS.hic_larvae_double_Rep1_dm6_BeS_1','hic.hic_larvae_double_Rep2_dm6_BeS.hic_larvae_double_Rep2_dm6_BeS_1','hic.hic_larvae_double_Rep3_dm6_BeS.hic_larvae_double_Rep3_dm6_BeS_1'),
	hic_larvae_DPRE1_merge_dm6_BeS          = c('hic.hic_larvae_DPRE1_Rep1_dm6_BeS.hic_larvae_DPRE1_Rep1_dm6_BeS_1','hic.hic_larvae_DPRE1_Rep2_dm6_BeS.hic_larvae_DPRE1_Rep2_dm6_BeS_1'),
	hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS = c('hic.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS_1','hic.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS_1')
	)

obsData <- list(
	hic_larvae_DWT_merge_dm6_BeS            = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1')
	)

print(refData)
print(obsData)
refTracks  <- refData
obsTracks <- obsData
print(refTracks)
print(obsTracks)

registerDoParallel(cores=12)

chr <- 'chr2L'
targetRegion2D <- g2d(gintervals(chr,0,gintervals.all()[gintervals.all()$chrom == chr,'end']))
print(targetRegion2D)

step <- 1e6
expand <- step/3
k    <- 250
kexp <- k
 
registerDoParallel(cores=12)

randMinTracks <- gtrack.ls('randMin')
print(randMinTracks)

chr   <- 'chr2L'
start <- 10e6
end   <- 20e6

chrIntrv <- gintervals(chr,start,end)
chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step))
print(chrIntrv)

#print("Computing differential maps!")

for(obsTrack in randMinTracks)
{
    print(obsTrack)
    for(expTrack in randMinTracks)
    {
        print(expTrack)

	if(obsTrack == expTrack){next}	

	hicRoot <- "hic.diffMaps."
 		
	set1 <- gsub('_randMin','',unlist(strsplit(obsTrack,'\\.'))[3])
	set2 <- gsub('_randMin','',unlist(strsplit(expTrack,'\\.'))[3])

	trackName <- paste0(set1,'_vs_',set2,'_',chr,'_score_k',k,'_kexp',kexp,"_",step/1e6,'Mb')
	print(trackName)
	outTrack <- paste0(hicRoot,trackName)		
	print(outTrack)
	hicName <- paste0(set1,'_vs_',set2)
	work_dir <- paste0("./_tmp_",trackName,"/")

	if(!dir.exists(work_dir)){dir.create(work_dir, mode="7777", recursive=TRUE);}

	if(gtrack.exists(outTrack))
	{
	    print(paste0(outTrack," exists"))
	    next
	}

	foreach(i=1:nrow(chrIter)) %dopar% {
	int1 <- chrIter[i,]
	int2 <- expand2D(int1,expand)[,1:6]
			
	outFile <- paste0(work_dir,chr,'_',i,'.scores')
	if(file.exists(outFile))
	{
	    print('done...')
   	    return(i)
	}
 			
	print('scoring...')
	scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2, k=k, k_exp=kexp)
	data <- scores$points
	write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)
 			
	return(i)
    }
	
    files <- list.files(work_dir, full.names=T,pattern='scores')
 		
    ### Single Import
    gtrack.2d.import(outTrack, paste("Score track for", hicName,' - ',chr), files)
    ### Add Int AND Reverse INT
    # gtrack.2d.import_contacts(outTrack, paste("import diffMaps scores for ", trackName), files)

    ### Remove temp files....
    for (f in files){file.remove(f)}
    }
}
quit()
