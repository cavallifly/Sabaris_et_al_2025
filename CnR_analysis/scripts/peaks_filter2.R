getwd()
setwd("N:/Commun-Cavalli/Gonzalo/2024_GAF_project/2024_Primary_analysis/macs3")

E_GAF_narrowPeaks <- read.delim("E_GAF_R1R2merged_macs3_q0.001_peaks.narrowPeak", header = FALSE, sep="\t")
E_GAF_narrowPeaks
nrow(E_GAF_narrowPeaks)
summary(E_GAF_narrowPeaks)

WD_GAF_narrowPeaks <- read.delim("WD_GAF_R1R2merged_macs3_q0.001_peaks.narrowPeak", header = FALSE, sep="\t")
nrow(WD_GAF_narrowPeaks)

#filter peaks by Interger Score > 100 (column #5) AND fold-change at peak summit >2 (column #7)
#5th: integer score for display. It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue)
E_GAF_narrowPeaks_s100_FC2 <- E_GAF_narrowPeaks[ which(E_GAF_narrowPeaks$V5 >100 & E_GAF_narrowPeaks$V7 > 2 ) , ]
nrow(E_GAF_narrowPeaks_s100_FC2)
write.table(E_GAF_narrowPeaks_s100_FC2, file="E_GAF_narrowPeaks_s100_FC2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2 <- WD_GAF_narrowPeaks[ which(WD_GAF_narrowPeaks$V5 >100 & WD_GAF_narrowPeaks$V7 > 2 ) , ]
nrow(WD_GAF_narrowPeaks_s100_FC2)
write.table(WD_GAF_narrowPeaks_s100_FC2, file="WD_GAF_narrowPeaks_s100_FC2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#quartiles
summary(E_GAF_narrowPeaks_s100_FC2)
E_GAF_narrowPeaks_s100_FC2$V11<- cut(E_GAF_narrowPeaks_s100_FC2$V5,quantile(E_GAF_narrowPeaks_s100_FC2$V5),include.lowest=TRUE,labels=FALSE)
E_GAF_narrowPeaks_s100_FC2
write.table(E_GAF_narrowPeaks_s100_FC2, file="E_GAF_narrowPeaks_s100_FC2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


E_GAF_narrowPeaks_s100_FC2_q1<- E_GAF_narrowPeaks_s100_FC2[ which(E_GAF_narrowPeaks_s100_FC2$V11 ==1) , ]
nrow(E_GAF_narrowPeaks_s100_FC2_q1)
write.table(E_GAF_narrowPeaks_s100_FC2_q1, file="E_GAF_narrowPeaks_s100_FC2_q1.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

E_GAF_narrowPeaks_s100_FC2_q2<- E_GAF_narrowPeaks_s100_FC2[ which(E_GAF_narrowPeaks_s100_FC2$V11 ==2) , ]
nrow(E_GAF_narrowPeaks_s100_FC2_q2)
write.table(E_GAF_narrowPeaks_s100_FC2_q2, file="E_GAF_narrowPeaks_s100_FC2_q2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

E_GAF_narrowPeaks_s100_FC2_q3<- E_GAF_narrowPeaks_s100_FC2[ which(E_GAF_narrowPeaks_s100_FC2$V11 ==3) , ]
nrow(E_GAF_narrowPeaks_s100_FC2_q3)
write.table(E_GAF_narrowPeaks_s100_FC2_q3, file="E_GAF_narrowPeaks_s100_FC2_q3.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

E_GAF_narrowPeaks_s100_FC2_q4<- E_GAF_narrowPeaks_s100_FC2[ which(E_GAF_narrowPeaks_s100_FC2$V11 ==4) , ]
nrow(E_GAF_narrowPeaks_s100_FC2_q4)
write.table(E_GAF_narrowPeaks_s100_FC2_q4, file="E_GAF_narrowPeaks_s100_FC2_q4.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2$V11<- cut(WD_GAF_narrowPeaks_s100_FC2$V5,quantile(WD_GAF_narrowPeaks_s100_FC2$V5),include.lowest=TRUE,labels=FALSE)
WD_GAF_narrowPeaks_s100_FC2
write.table(WD_GAF_narrowPeaks_s100_FC2, file="WD_GAF_narrowPeaks_s100_FC2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2_q1<- WD_GAF_narrowPeaks_s100_FC2[ which(WD_GAF_narrowPeaks_s100_FC2$V11 ==1) , ]
nrow(WD_GAF_narrowPeaks_s100_FC2_q1)
write.table(WD_GAF_narrowPeaks_s100_FC2_q1, file="WD_GAF_narrowPeaks_s100_FC2_q1.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2_q2<- WD_GAF_narrowPeaks_s100_FC2[ which(WD_GAF_narrowPeaks_s100_FC2$V11 ==2) , ]
nrow(WD_GAF_narrowPeaks_s100_FC2_q2)
write.table(WD_GAF_narrowPeaks_s100_FC2_q2, file="WD_GAF_narrowPeaks_s100_FC2_q2.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2_q3<- WD_GAF_narrowPeaks_s100_FC2[ which(WD_GAF_narrowPeaks_s100_FC2$V11 ==3) , ]
nrow(WD_GAF_narrowPeaks_s100_FC2_q3)
write.table(WD_GAF_narrowPeaks_s100_FC2_q3, file="WD_GAF_narrowPeaks_s100_FC2_q3.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WD_GAF_narrowPeaks_s100_FC2_q4<- WD_GAF_narrowPeaks_s100_FC2[ which(WD_GAF_narrowPeaks_s100_FC2$V11 ==4) , ]
nrow(WD_GAF_narrowPeaks_s100_FC2_q4)
write.table(WD_GAF_narrowPeaks_s100_FC2_q4, file="WD_GAF_narrowPeaks_s100_FC2_q4.narrowPeak", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
