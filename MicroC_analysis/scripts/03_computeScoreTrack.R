require(misha)
require(shaman)

require(doParallel)

args = commandArgs(trailingOnly=TRUE)

##########################################
### misha working DB
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('./scripts/auxFunctions.R')
options(scipen=20,gmax.data.size=0.6e8,shaman.sge_support=1)

chr <- 'chr2L'
#start <- 10e6
#end   <- 20e6

step   <- 1000000 # Size of the chunk to divide the calculation
expand <- step/3  # Size of skin around the main chunk to look for expected counts v2
k      <- 250     # Number of neighbours to look for
k_exp  <- 2*k 

samples = c("larvae_DWT","larvae_double","larvae_DPRE1Up")
print(samples)	

for(sample in samples)
{
    print(sample)	
    if(sample == "larvae_DWT")
    {
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep1_dm6_BeS/hic_larvae_DWT_Rep1_dm6_BeS_1.track:
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep2_dm6_BeS/hic_larvae_DWT_Rep2_dm6_BeS_1.track:
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DWT_Rep3_dm6_BeS/hic_larvae_DWT_Rep3_dm6_BeS_1.track:
	nreplicates = 3
	author      = "BeS"
    }

    if(sample == "larvae_double")
    {
    	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_double_Rep1_dm6_BeS/hic_larvae_double_Rep1_dm6_BeS_1.track
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_double_Rep2_dm6_BeS/hic_larvae_double_Rep2_dm6_BeS_1.track
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_double_Rep3_dm6_BeS/hic_larvae_double_Rep3_dm6_BeS_1.track	
	nreplicates = 3
	author      = "BeS"
    }

    if(sample == "larvae_DPRE1")
    {
    	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DPRE1_Rep1_dm6_BeS/hic_larvae_DPRE1_Rep1_dm6_BeS_1.track
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DPRE1_Rep2_dm6_BeS/hic_larvae_DPRE1_Rep2_dm6_BeS_1.track
	nreplicates = 2
	author      = "BeS"
    }

    if(sample == "larvae_DPRE1Up")
    {
    	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DPRE1Up_Rep1_dm6_BeS/hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS_1.track
	#/zdata/data/mishaDB/trackdb/dm6/tracks/hic/hic_larvae_DPRE1Up_Rep2_dm6_BeS/hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS_1.track
	nreplicates = 2
	author      = "BeS"
    }    

    for(chr in c(chr))
    #for(chr in gintervals.all()$chrom)    
    {
        print(chr)
	start <- gintervals.all()[gintervals.all()$chrom == chr,]$start
	end   <- gintervals.all()[gintervals.all()$chrom == chr,]$end	

	chrIntrv <- gintervals(chr,start,end)	
	chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step)) #,band=c(-10e10,0))
	print(chr);print(chrIntrv$start);print(chrIntrv$end)
	print(nrow(chrIter))

	### SCORE MAPS
	if(nreplicates == 2)
    	{
	    obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"))
	    if (sample == "larvae_DPRE1Up")
	    {
	        obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6DPRE1Up_",author,".hic_",sample,"_Rep1_dm6DPRE1Up_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6DPRE1Up_",author,".hic_",sample,"_Rep2_dm6DPRE1Up_",author,"_1"))
	    }
	}
	if(nreplicates == 3)
	{        
            obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep3_dm6_",author,".hic_",sample,"_Rep3_dm6_",author,"_1"))
	    if(sample == "LE_WT")
	    {
	        obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_YO.hic_",sample,"_Rep1_dm6_YO_1"),paste0("hic.hic_",sample,"_Rep2_dm6_YO.hic_",sample,"_Rep2_dm6_YO_1"),paste0("hic.hic_",sample,"_Rep3_dm6_VL.hic_",sample,"_Rep3_dm6_VL_1"))
	    }
        }
        if(nreplicates == 4)
        {
            obsTrack <- c(paste0("hic.hic_",sample,"_Rep1_dm6_",author,".hic_",sample,"_Rep1_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2_dm6_",author,".hic_",sample,"_Rep2_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep1res_dm6_",author,".hic_",sample,"_Rep1res_dm6_",author,"_1"),paste0("hic.hic_",sample,"_Rep2res_dm6_",author,".hic_",sample,"_Rep2res_dm6_",author,"_1"))
        }

        setName  <- paste0("hicScores_",sample,"_merge_dm6_",author)
        if(sample == "LE_WT")
	{
	    if(author == "YOVL")
	    {
                setName  <- paste0("hicScores_",sample,"_merge_dm6_YOVL")
	    }
        }

	expTrack <- paste0(obsTrack,"_shuffle_500Small_1000High")
	print(obsTrack)
	print(setName)
	print(expTrack)

	hicName <- paste0("hic")
	
	trackName <- paste0(setName,"_",chr,"_",chrIntrv$start,"_",chrIntrv$end,"bp_score_k",k,"_kexp",k_exp,"_step",step/1e6,"Mb")
	outTrack <- paste0(hicName,".",setName,".",trackName)
	print(outTrack)

	work_dir <- paste0("./_tmp_",trackName,"/")
	if(!dir.exists(work_dir))
	{
	    dir.create(work_dir, mode="7777", recursive=TRUE)
        } #else {
	#    next
	#}

	if(gtrack.exists(outTrack))
	{
            print(paste0(outTrack," exists!"))
            next
        }

        registerDoParallel(cores=5)

	foreach(i=1:nrow(chrIter)) %dopar% {
     	    int1 <- chrIter[i,]	
	    int2 <- expand2D(int1,expand)[,1:6]
	    print(paste0("Interval to score"))
	    print(int1)	
	    print(paste0("Interval + skin"))	
	    print(int2)

	    outFile <- paste0(work_dir,chr,"_",i,".scores")
	    if(file.exists(outFile))
	    {
	        print("done...")
		return(i)
	    }

	    genDist <- abs(int1$start2-int1$start1)
	    if(genDist <= 500e3)
	    {
	        print(paste0("Closer than ",500e3,"bp. Leave it for later!"))	
	    	return(i)
	    }

	    print(paste0("scoring portion ",i,"..."))
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k,k_exp=k_exp)
	    data <- scores$points
	    write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	    return(i)
	}

	registerDoParallel(cores=1)	
	
	foreach(i=1:nrow(chrIter)) %dopar% {
    	    int1 <- chrIter[i,]	
	    int2 <- expand2D(int1,expand)[,1:6]
	    print(paste0("Interval to score"))
	    print(int1)	
	    print(paste0("Interval + skin"))	
	    print(int2)

	    outFile <- paste0(work_dir,chr,"_",i,".scores")
	    if(file.exists(outFile))
	    {
	        print("done...")
	    	return(i)
	    }

	    genDist <- abs(int1$start2-int1$start1)
	    if(genDist > 500e3)
	    {
		print(paste0("Further than ",500e3,"bp. Already done!"))
		return(i)
	    }

	    print(paste0("scoring portion ",i,"..."))
	    scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k,k_exp=k_exp)
	    data <- scores$points
	    write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	    return(i)
	}

	files <- list.files(work_dir, full.names=T,pattern="scores")
	
	trackFolder <- paste0(dbDir,"tracks/",gsub("\\.","/",hicName),"/",setName,"/")
	if(!dir.exists(trackFolder)){dir.create(trackFolder, mode="7777", recursive=TRUE);}

	### Single Import
	gtrack.2d.import(outTrack, paste("Score track for", hicName," - ",chr), files)
	### Add Int AND Reverse INT
	# gtrack.2d.import_contacts(outTrack, paste("import diffMaps scores for ", trackName), files)

	### Remove temp files....
	for (f in files){file.remove(f)}
    }
}
