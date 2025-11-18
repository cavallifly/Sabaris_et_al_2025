#!/bin/env Rscript

library(misha)
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()
source('./scripts_clean/auxFunctions.R')

obsData <- list(
	microc_ED_PH18_merge_dm6_BeS = c('microC.mC_ED_PH18_Rep1_dm6_BeS.mC_ED_PH18_Rep1_dm6_BeS', 'microC.mC_ED_PH18_Rep1res_dm6_BeS.mC_ED_PH18_Rep1res_dm6_BeS', 'microC.mC_ED_PH18_Rep2_dm6_BeS.mC_ED_PH18_Rep2_dm6_BeS', 'microC.mC_ED_PH18_Rep2res_dm6_BeS.mC_ED_PH18_Rep2res_dm6_BeS'),
	microc_ED_PH29_merge_dm6_BeS = c('microC.mC_ED_PH29_Rep1_dm6_BeS.mC_ED_PH29_Rep1_dm6_BeS', 'microC.mC_ED_PH29_Rep1res_dm6_BeS.mC_ED_PH29_Rep1res_dm6_BeS', 'microC.mC_ED_PH29_Rep2_dm6_BeS.mC_ED_PH29_Rep2_dm6_BeS', 'microC.mC_ED_PH29_Rep2res_dm6_BeS.mC_ED_PH29_Rep2res_dm6_BeS')
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
