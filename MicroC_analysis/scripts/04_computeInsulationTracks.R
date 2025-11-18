#!/bin/env Rscript

library(misha)
mDBloc <-  './mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()
source('./scripts/auxFunctions.R')

obsData <- list(
	hic_larvae_DWT_merge_dm6_BeS            = c('hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1','hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1','hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1'),
	hic_larvae_double_merge_dm6_BeS         = c('hic.hic_larvae_double_Rep1_dm6_BeS.hic_larvae_double_Rep1_dm6_BeS_1','hic.hic_larvae_double_Rep2_dm6_BeS.hic_larvae_double_Rep2_dm6_BeS_1','hic.hic_larvae_double_Rep3_dm6_BeS.hic_larvae_double_Rep3_dm6_BeS_1'),
	hic_larvae_DPRE1_merge_dm6_BeS          = c('hic.hic_larvae_DPRE1_Rep1_dm6_BeS.hic_larvae_DPRE1_Rep1_dm6_BeS_1','hic.hic_larvae_DPRE1_Rep2_dm6_BeS.hic_larvae_DPRE1_Rep2_dm6_BeS_1'),
	hic_larvae_DPRE1Up_merge_dm6DPRE1Up_BeS = c('hic.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep1_dm6DPRE1Up_BeS_1','hic.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS.hic_larvae_DPRE1Up_Rep2_dm6DPRE1Up_BeS_1')
	)

for(k in seq(100, 300, by=50))
{
    insScale <- k * 1e+3
    insResolution <- 2e+3
	
    for(set in names(obsData))
    {
	insTrack     <- obsData[[set]]
	insTrackName <- paste0('insulation.INS_',set,'_w',insScale/1000,'kb','_r',insResolution/1000,'kb')
	print(insTrackName)
	print(gtrack.exists(insTrackName))
	if(gtrack.exists(insTrackName)){print(paste0("Track ",insTrackName," exists!")); next}
	print(insTrack)
	gtrack.2d.gen_insu_track(insTrack, insScale, insResolution, min_diag_d=1000, insTrackName, description="")
    }
}
