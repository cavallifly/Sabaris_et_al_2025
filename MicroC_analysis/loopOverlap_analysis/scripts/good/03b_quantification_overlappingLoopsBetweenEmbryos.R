# Script inspired by https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova

library(tidyverse)
library(ggpubr)
library(rstatix)
library(dunn.test)
library(dplyr)
library(FSA)

inFiles = list.files("./",pattern=".*BeS.*_range4.*_res250_centralRegion_.*pixels_data.txt.*")
print(inFiles)

subsetting <- "FALSE"

colors=c(rgb( 95, 158, 160, 1.0, names = NULL, maxColorValue = 255), rgb( 106, 90, 205, 1.0, names = NULL, maxColorValue = 255), rgb( 128, 128, 128, 1.0, names = NULL, maxColorValue = 255))

for(inFile in inFiles)
{
    print(inFile, quote=F)
    allData <- read.table(inFile, header=F)
    colnames(allData) <- c("condition","values")    
    print(head(allData), quote=F)
    name <- gsub(".txt","",inFile)
    name <- gsub(".tsv","",name)        

    data <- allData

    levels <- unique(data$condition)
    print(levels)
    data$condition = factor(data$condition, levels=levels)
    print(data)

    print(paste0("### Get summary statistics ###"), quote=F)
    dataStats <- data %>%
	  group_by(condition) %>%
	  get_summary_stats(values, type = "mean_sd")
    print(as.data.frame(dataStats), quote=F)

    summ <- data %>%
           group_by(condition) %>%
	   summarize(n = n(), condition=condition, score = ((1-0.05)*min(data$values)))
    summ <- unique(summ)
    print(as.data.frame(summ), quote=F)   

    bxp <- ggviolin(data, x="condition", y="values", color="black", fill="condition", alpha=.50, scale = "width", trim = TRUE) + scale_fill_manual(values=colors) +
           geom_boxplot(position=position_dodge(1), width=0.1, color="black", outlier.shape = NA) +
           geom_point(position = position_jitter(width=0.4), size = 0.5, color="black", alpha=1.0) +	   	   
	   labs(y="log2(observed/expected)", x = "") +
	   geom_text(data=summ, aes(x=condition, y=score, label = paste0("n = ",n)), size=4., color="black") +    
	   theme_classic() +
   	   theme(axis.text.x = element_text(angle = 60, hjust=1, size=14), axis.text.y = element_text(size=14), legend.position = "none",
	   axis.title.y = element_text(size=16))

    res.aov <- kruskal_test(values ~ condition, data = data)
    print(res.aov, quote=F)
    significantInteraction = res.aov[res.aov$p < 0.05,]
    if(nrow(significantInteraction) > 0)
    {
        print(paste0("There was a statistically significant interaction between locis"), quote=F)
        print(significantInteraction, quote=F)
    } else {
        print(paste0("There was no statistically significant interaction between loci"), quote=F)
    }

    print("", quote=F)	
    print(paste0("### Post-hoct tests ###"), quote=F)
    pairwiseTestsCondition <- dunn_test(values ~ condition, data = data, p.adjust.method = "BH")

    print(paste0("### Visualization: violin plots with p-values"), quote=F)

    print(as.data.frame(pairwiseTestsCondition), quote=F)
    outFileText <- paste0("violinPlot_",name,"_with_stats.tsv")
    print(outFileText)
    write.table(as.data.frame(pairwiseTestsCondition), file=outFileText,  row.names=F, sep="\t")    
    pairwiseTestsCondition <- pairwiseTestsCondition %>% 
    			      			       add_xy_position(x = "condition")

    finalBxp <- bxp +
	   stat_pvalue_manual(data=pairwiseTestsCondition, label="p.adj.formatted", label.size = 4, bracket.size = 0.2, tip.length = 0, ) +
	   labs(subtitle = get_test_label(res.aov, detailed = TRUE), caption = get_pwc_label(pairwiseTestsCondition)) +
	   theme(legend.position = "none",
               plot.title = element_text(size = 18, color = "black", face= "bold"),
               plot.subtitle = element_text(size = 16, face = "bold.italic"),
	       plot.caption = element_text(size = 16, face = "bold.italic")	       
               )

    outFile <- paste0("violinPlot_",name,"_with_stats.pdf")
    pdf(outFile)
    print(finalBxp)
    dev.off()	
    print("", quote=F)

} # Close cycle over inFile