library(ComplexUpset)
library(tidyr)
library(UpSetR)

args = commandArgs(trailingOnly=TRUE)

slop <- args[[1]]

inFiles = list.files("./",pattern=paste0(".*_renamedLoops_ext_",slop,"bp.bedpe",".*"))
#inFiles = inFiles[-grep("all",inFiles)]
print(inFiles)
#quit()

sets <- list()
for(inFile in inFiles)
{
    name = gsub("Loops_","",gsub(paste0("_renamedLoops_ext_",slop,"bp.bedpe"),paste0("_ext_",as.integer(slop)/1000,"kb"),inFile))
    print(paste0(inFile," ",name))

    loops <- read.table(inFile, header=F)
    colnames(loops) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "loop")
    print(head(loops))

    sets[[name]] <- loops$loop
    print(sets)
}

outFile = paste0("upset_plot_loop_overlap.pdf")
p <- upset(fromList(sets), order.by="freq", nintersects=NA)
pdf(outFile)
print(p)
dev.off()
