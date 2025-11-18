#!/bin/env Rscript

library(misha)
library(shaman)

require(doParallel)
registerDoParallel(cores=4)

options(shaman.mc_support=1, gmax.data.size=2.4e8)

assembly <- 'dm6'
db       <- paste0('./mishaDB/trackdb/',assembly,'/')
gsetroot(db)
gdb.reload()

obsTracks <- gtrack.ls('hic.hic_')
print(obsTracks)

foreach(i=1:length(obsTracks)) %dopar%
{
    currentTrack <- obsTracks[i]
	
    if(length(grep('shuffle|randMin',currentTrack)) > 0){next}
    hicName <- unlist(strsplit(currentTrack,'\\.'))[3]
    #print(hicName)

    grdSmall <- 500e3
    grdHigh <- 1000e3
    outName <- paste0(currentTrack,'_shuffle_',grdSmall/1e3,'Small_',grdHigh/1e3,'High')
    print(paste0("outName ",outName))

    if(gtrack.exists(outName))
    {
        message(currentTrack,' already done...')
        return()
    }
    message('Working on: ',currentTrack)
	
    tmpDir="./tempData/"
    if(!dir.exists(tmpDir)){dir.create(tmpDir, mode="7777", recursive=TRUE);}
    currentWD <- paste0(tmpDir,hicName,'_shuffle_',grdSmall/1e3,'Small_',grdHigh/1e3,'High')
    print(currentWD)
    message('Temp data in: ',currentWD)

    if(!dir.exists(currentWD)){dir.create(currentWD)}
	
    shaman_shuffle_hic_track(track_db=db, obs_track_nm=currentTrack, exp_track_nm=outName, work_dir=currentWD, max_jobs=2, grid_small=grdSmall, grid_high=grdHigh)
	
    ### if import fails... sed -i '1,5s/-e //' *.uniq
    files <- list.files(currentWD, full.names=T,pattern='uniq')
    gtrack.2d.import(outName, paste("shuffled 2d track for", currentTrack), files)
    return()
}
