library(scales) # to access break formatting functions
library(ggplot2)
library(dplyr)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

h = as.numeric(args[1])

data <- read.table('_tmp.txt')
colnames(data) <- c("condition","score")
print(head(data))

data$condition = factor(data$condition, levels=c("Control","PH_KD"))
colors=c(rgb( 0,   0,    225, 1.0, names = NULL, maxColorValue = 255), rgb( 83, 83, 83, 1.0, names = NULL, maxColorValue = 255))
#colors=c("#CF6263","#F9E86F")
print(head(data))
#PH18 PH29 PHD11

summ <- data %>%
  group_by(across(all_of(c("condition")))) %>%
  summarize(n = n(), score = h)
print(summ)

pdf('_tmp1.pdf')

compare_means(score ~ condition,  data = data, method = "wilcox.test", group_by="condition")
print(compare_means)

my_comparisons <- list( c("Control", "PH_KD") )
plot <- ggplot(data, aes(x=condition, y=score, color = NA, fill=condition), alpha=0.7) + geom_violin(position=position_dodge(1), trim=T, scale="width") + geom_boxplot(position=position_dodge(1), width=0.1, color="black", outlier.shape = NA) + scale_color_manual(values=colors) + scale_fill_manual(values=colors) + theme(panel.background = element_rect(fill = NA), panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), axis.title = element_text(face="bold",size=24), axis.text.x = element_text(angle=60, hjust=1), axis.text = element_text(face="bold",size=11), axis.ticks.x = element_blank(), legend.position = "top", legend.title =  element_blank()) + labs(x="",y="log2(observed/expected)") + geom_text(aes(label=n), color="black", data = summ, position=position_dodge(1)) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + geom_text(aes(label=n), color="black", data = summ, position=position_dodge(1)) + guides(color="none")

#stat_compare_means(aes(group = condition), method = "wilcox.test") + geom_text(aes(label=n), color="black", data = summ, position=position_dodge(1)) + guides(color="none")
#panel.grid.major.y = element_line(colour = "grey10"),
print(plot)

#+ stat_compare_means(clabel = "p.format")
dev.off()

#+ stat_summary(fun.data = give.n, geom = "text")


#+ guides(file="none")
#
#+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits= c(1000,3000000)) + annotation_logticks(sides="l")

#geom_text(data=Summary.data, aes(x = V1, y = V2, label=n),color="red", fontface =2, size = 5)
#theme(axis.title = element_text(face="bold",size=24), axis.text = element_text(face="bold",size=11), axis.text.x = element_text(angle=60, hjust=1)) + labs(x="",y="Genomic Distance (bp)",fill="Distances")


#dev.off()
