library(ComplexUpset)
library(tidyr)
library(UpSetR)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

slop <- args[[1]]

inFiles = list.files("./",pattern=paste0(".*_renamedLoops_ext_",slop,"bp.bedpe",".*"))
#inFiles = inFiles[-grep("all",inFiles)]
print(inFiles)
#quit()

sets <- list()
for(inFile in inFiles)
{
    name = gsub("Loops_Embryo_WT_","",gsub(paste0("_renamedLoops_ext_",slop,"bp.bedpe"),"",inFile))
    name = gsub("BeS","Sabaris2026",name)
    print(paste0(inFile," ",name))

    loops <- read.table(inFile, header=F)
    colnames(loops) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "loop")
    print(head(loops))

    sets[[name]] <- loops$loop
    print(sets)
}
print(names(sets))
#print(names(sets)[1])
#print(names(sets)[2])
print(colnames(fromList(sets)))

outFile = paste0("upset_plot_loop_overlap.pdf")
upset_plot <- upset(fromList(sets), nintersects=NA, text.scale = 1.5, order.by = "freq")

p <- upset_plot #+
     #theme(axis.text.x = element_text(size = 12))
pdf(outFile)
print(p)
dev.off()
