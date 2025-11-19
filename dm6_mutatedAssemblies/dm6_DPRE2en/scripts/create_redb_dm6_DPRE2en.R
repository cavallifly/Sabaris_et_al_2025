#!/bin/env Rscript

rm(list=ls())
rLibLoc <- '/home/michael.szalay/anaconda3/envs/misha/lib/R/library/'
#######################################################
require('misha', lib.loc=rLibLoc, quietly=TRUE)

setwd("/home/minHiC2_sourceFiles/auxScripts/")
source("TG3C/imp3c.r")

### dbDir ### trackName ### runConfig
dbDir <- '/zdata/data/mishaDB/trackdb/dm6_DPRE2en/'

gsetroot(dbDir)
gdb.reload()

#### redb track
restrictionSeq <- 'GATC'

checkSp <- gsub('/','',dbDir)
#print(checkSp)
checkSp <- unlist(strsplit(checkSp,'trackdb'))[2]
print(checkSp)
#quit()

if(length(gtrack.ls(restrictionSeq)) < 1)
{
    # First one as to create the mapability track
    gtrack.create_mapab_track(paste0("/zdata/data/2024_01_08_HiC_Sandrine2/00_build_the_dm6_DPRE2en_modified_genome/mapa_",checkSp,".config"))
    # Second one as to creat the restriction enzyme database track
    gtrack.create_redb_tracks(restrictionSeq,paste0('/zdata/data/2024_01_08_HiC_Sandrine2/00_build_the_dm6_DPRE2en_modified_genome/redb_',checkSp,'.conf'),verbose=TRUE)
}


gsetroot(dbDir)
gdb.reload()

