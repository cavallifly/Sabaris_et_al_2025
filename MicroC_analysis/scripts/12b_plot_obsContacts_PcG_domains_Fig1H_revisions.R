require(reshape2) # cat tabulate to matrix
library(ggplot2)
library(dplyr)
library(ggpubr)


args = commandArgs(trailingOnly=TRUE)
source('./scripts_clean/auxFunctions.R')
options("scipen"=999)
colors <- rainbow(10)

obsData <- list(
	larvae_double  = c("hic.hic_larvae_double_Rep1_dm6_BeS.hic_larvae_double_Rep1_dm6_BeS_1",
			   "hic.hic_larvae_double_Rep2_dm6_BeS.hic_larvae_double_Rep2_dm6_BeS_1",
			   "hic.hic_larvae_double_Rep3_dm6_BeS.hic_larvae_double_Rep3_dm6_BeS_1"),
	LD_gypsy2Enh   = c("hic.hic_LD_gypsy2Enh_Rep1_dm6_BeS.hic_LD_gypsy2Enh_Rep1_dm6_BeS_1",
   		           "hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1",
   		           "hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1")
		       )
samples = names(obsData)
print(paste0("Samples ",samples))

refSample <- "larvae_DWT"
refChrom  <- "chr2L"

outfile <- paste0("obsContacts_within_domains_vs_",refSample,"_in_",refChrom,"_revisions.tab")
df <- read.table(outfile,header=T)
print(head(df))
toMatch <- c("chr2L")
df <- df[grep(paste(toMatch,collapse="|"),df$chrom1),]
df <- df[-grep("larvae_DWT",df$sample),]
print(head(df))

summ <- df              %>%
    group_by(sample)  %>%
    #reframe(n = n(), sample=unique(sample), y=min(df$Log2FC_Contacts_larvae_DWT-0.025))
    reframe(n = n(), sample=unique(sample), y=min(df$Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_larvae_DWT-0.025))    

levels <- samples
my_comparisons <- list(c("larvae_DWT","larvae_double"),c("larvae_DWT","LD_gypsy2Enh"))

df$sample<- factor(df$sample, levels = levels)

print(head(df))
#quit()

#statAnalysis <- compare_means(Log2FC_Contacts_larvae_DWT ~ sample, data = df, method = "wilcox.test")
statAnalysis <- compare_means(Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_larvae_DWT ~ sample, data = df, method = "wilcox.test")
write.table(statAnalysis, file = "statAnalysis_normObsContacts_between_PcG_domains.tab", sep="\t", row.names=FALSE, quote=FALSE)

#p <- ggplot(df, aes(x=sample, y=Log2FC_Contacts_larvae_DWT, fill=sample, alpha=0.2)) +
p <- ggplot(df, aes(x=sample, y=Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_larvae_DWT, fill=sample, alpha=0.2)) +
     geom_violin(position = position_dodge(1), trim=T, scale = "width") +
     geom_boxplot(position = position_dodge(1), width=0.1, fill="white", color="black", outlier.shape = NA) +
     #ylim(c(min(df$Log2FC_Contacts_larvae_DWT)-0.05,max(df$Log2FC_Contacts_larvae_DWT))) +
     ylim(c(min(df$Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_larvae_DWT)-0.05,max(df$Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_larvae_DWT))) +     
     theme(panel.background = element_rect(fill = NA),
           panel.grid.major.x = element_blank(),
           panel.grid.major.y = element_blank(),
           axis.title = element_text(size=12),	   
           axis.text = element_text(face="bold",size=11),
           axis.text.x = element_text(angle=60, hjust=1),
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "none") +
     labs(x="",y="Log2(Cov. norm. Mut_intraTAD/Cov. norm. larvae_DWT) intraTAD_counts)",fill="") +	   
     scale_color_manual(values=colors) +
     scale_fill_manual(values=colors) +
     #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p", vjust = 0.1) +
     geom_point(data=df[df$domain1 == "dac",],color="red",alpha=1) +
     geom_jitter(data=df[df$domain1 != "dac",],color="black",width=0.2) +
     geom_text(data=summ, aes(x=sample, y=y, label=paste0("n=",n))) +
     guides(fill="none")

pdf(paste0("obsContacts_within_domains_vs_",refSample,"_in_",refChrom,"_Fig1H_revisions.pdf"))
print(p)
dev.off()
