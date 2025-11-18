require(misha)
require(shaman)

require(doParallel)

### misha working DB
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source("/zdata/data/auxFunctions/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

obsData <- list(
	ED_10kb              = "ED_10kb",
	ED_10kb_new_internal = "ED_10kb_new_internal",
	ED_10kb_new	     = "ED_10kb_new",
	ED_10kb_withGaps     = "ED_10kb_withGaps",
	ED_10kb_withLoops    = "ED_10kb_withLoops",
	EDph_10kb	     = "EDph_10kb",
	EDph_10kb_new	     = "EDph_10kb_new",
	EDwt_10kb_new	     = "EDwt_10kb_new"
)


wDir = "domains."



for( set in obsData)
{
    tssCoordinates <- gintervals.load(paste0(wDir,obsData[[set]]))
    print(tssCoordinates)
    #rownames(tssCoordinates) <- tssCoordinates$geneName
}