library(tidyverse)
library(rstatix)
library(ggpubr)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(hrbrthemes)
library(viridis)

# setwd("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/FIMO_MEME/")
setwd("D:/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/FIMO_MEME/ED_peaks")

### ED peaks motifs
ED_GAF_Q1_motifs <- read.table('MA02052/ED_q1/fimo.tsv', header=T)
ED_GAF_Q2_motifs <- read.table('MA02052/ED_q2/fimo.tsv', header=T)
ED_GAF_Q3_motifs <- read.table('MA02052/ED_q3/fimo.tsv', header=T)
ED_GAF_Q4_motifs <- read.table('MA02052/ED_q4/fimo.tsv', header=T)

#filter table by q-value <=0.05
ED_GAF_Q1_motifs <- ED_GAF_Q1_motifs[ED_GAF_Q1_motifs$q.value <= 0.05,]
ED_GAF_Q2_motifs <- ED_GAF_Q2_motifs[ED_GAF_Q2_motifs$q.value <= 0.05,]
ED_GAF_Q3_motifs <- ED_GAF_Q3_motifs[ED_GAF_Q3_motifs$q.value <= 0.05,]
ED_GAF_Q4_motifs <- ED_GAF_Q4_motifs[ED_GAF_Q4_motifs$q.value <= 0.05,]


ED_GAF_Q1_motifs <- ED_GAF_Q1_motifs[order(ED_GAF_Q1_motifs$sequence_name),]
ED_GAF_Q2_motifs <- ED_GAF_Q2_motifs[order(ED_GAF_Q2_motifs$sequence_name),]
ED_GAF_Q3_motifs <- ED_GAF_Q3_motifs[order(ED_GAF_Q3_motifs$sequence_name),]
ED_GAF_Q4_motifs <- ED_GAF_Q4_motifs[order(ED_GAF_Q4_motifs$sequence_name),]

# counting number of motifs
ED_GAF_Q1_motifs_count <- ED_GAF_Q1_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q2_motifs_count <- ED_GAF_Q2_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q3_motifs_count <- ED_GAF_Q3_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q4_motifs_count <- ED_GAF_Q4_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')


### random WD peaks
random_ED <- read.table('MA02052/ED_random/fimo.tsv' , header=T)
random_ED <- random_ED[random_ED$q.value <= 0.05,]
random_ED <- random_ED[order(random_ED$sequence_name),]
random_ED_count <- random_ED %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')


#stats score
#prepare table for stats
ED_random <- cbind(random_ED, class= "random")
ED_Q1 <- cbind(ED_GAF_Q1_motifs, class= "Q1")
ED_Q2 <- cbind(ED_GAF_Q2_motifs, class= "Q2")
ED_Q3 <- cbind(ED_GAF_Q3_motifs, class= "Q3")
ED_Q4 <- cbind(ED_GAF_Q4_motifs, class= "Q4")
ED_score_stats <- rbind(ED_random, ED_Q1, ED_Q2, ED_Q3, ED_Q4)


ED_score_stats %>% group_by(class) %>% get_summary_stats(score, type = "median_iqr")
 
# kruskal_test with Dunnets for multiple comparisons
ED_score_kruskal <- ED_score_stats %>% rstatix::kruskal_test (score ~ class)
ED_score_kruskal
#effect size
ED_score_stats %>% kruskal_effsize(score ~ class) #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# Pairwise comparisons using Dunn's test
pwc <- ED_score_stats %>% dunn_test(score ~ class, p.adjust.method = "BH") 
pwc
fwrite(pwc, file ="MA02052/ED_DunnTest_BHpadj_scores_q005.tsv",sep= "\t", quote = F)
#end stats score

#calculate the total number of motifs in each class and create a vector.
ED_total_motifs <- c(sum(random_ED_count$total_count), sum(ED_GAF_Q1_motifs_count$total_count), sum(ED_GAF_Q2_motifs_count$total_count), sum(ED_GAF_Q3_motifs_count$total_count), sum(ED_GAF_Q4_motifs_count$total_count))
ED_total_motifs
names(ED_total_motifs) <- c("random", "GAF_Q1", "GAF_Q2", "GAF_Q3", "GAF_Q4")

# Convert to data frame
df <- data.frame(
  Group = factor(names(ED_total_motifs), levels = names(ED_total_motifs)),
  Count = as.numeric(ED_total_motifs)
)

v <- ggplot(df, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.9, direction = -1) +
  labs(title = "Eye disc", y = "Total GAF motifs", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(ylim = c(0, 8500))
v
ggsave("MA02052/ED_Q1toQ4_loopAnchors_Total_n_motif_q0.05.pdf", plot = v, width = 7, height = 5)


#WD_random peaks with motifs
gr_random_ED_motifs <- GRanges(
  seqnames = random_ED$sequence_name,
  ranges = IRanges(start = random_ED$start,
                   end = random_ED$stop),
)

random_ED_peaks <- read.table('ED_random_bed_l967_n1610_seed24_NoChr.bed', header=F)
random_ED_peaks <- random_ED_peaks[order(random_ED_peaks$V1, random_ED_peaks$V2, random_ED_peaks$V3),]

gr_random_ED_peaks <- GRanges(
  seqnames = random_ED_peaks$V1,
  ranges = IRanges(start = random_ED_peaks$V2,
                   end = random_ED_peaks$V3),
)

overlap_random <- subsetByOverlaps(gr_random_ED_peaks, gr_random_ED_motifs)
mcols(gr_random_ED_peaks)$overlap <- countOverlaps(gr_random_ED_peaks, gr_random_ED_motifs)

#n of peaks with motifs
gr_random_ED_peaks[ gr_random_ED_peaks$overlap >= 1 ]
df_random_ED_peaks_overlapMotifs<- as.data.frame(gr_random_ED_peaks[ gr_random_ED_peaks$overlap >= 1 ])

# WD_Qn peaks with motifs
gr_Q1_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q1_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q1_motifs$start,
                   end = ED_GAF_Q1_motifs$stop),
)

gr_Q2_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q2_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q2_motifs$start,
                   end = ED_GAF_Q2_motifs$stop),
)

gr_Q3_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q3_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q3_motifs$start,
                   end = ED_GAF_Q3_motifs$stop),
)

gr_Q4_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q4_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q4_motifs$start,
                   end = ED_GAF_Q4_motifs$stop),
)

ED_Q1_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q1_NoChr.bed', header=F)
ED_Q1_peaks <- ED_Q1_peaks[order(ED_Q1_peaks$V1, ED_Q1_peaks$V2, ED_Q1_peaks$V3),]

ED_Q2_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q2_NoChr.bed', header=F)
ED_Q2_peaks <- ED_Q2_peaks[order(ED_Q2_peaks$V1, ED_Q2_peaks$V2, ED_Q2_peaks$V3),]

ED_Q3_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q3_NoChr.bed', header=F)
ED_Q3_peaks <- ED_Q3_peaks[order(ED_Q3_peaks$V1, ED_Q3_peaks$V2, ED_Q3_peaks$V3),]

ED_Q4_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q4_NoChr.bed', header=F)
ED_Q4_peaks <- ED_Q4_peaks[order(ED_Q4_peaks$V1, ED_Q4_peaks$V2, ED_Q4_peaks$V3),]

gr_ED_Q1_peaks <- GRanges( 
  seqnames = ED_Q1_peaks$V1,
  ranges = IRanges(start = ED_Q1_peaks$V2,
                   end = ED_Q1_peaks$V3),
)

gr_ED_Q2_peaks <- GRanges( 
  seqnames = ED_Q2_peaks$V1,
  ranges = IRanges(start = ED_Q2_peaks$V2,
                   end = ED_Q2_peaks$V3),
)

gr_ED_Q3_peaks <- GRanges( 
  seqnames = ED_Q3_peaks$V1,
  ranges = IRanges(start = ED_Q3_peaks$V2,
                   end = ED_Q3_peaks$V3),
)

gr_ED_Q4_peaks <- GRanges( 
  seqnames = ED_Q4_peaks$V1,
  ranges = IRanges(start = ED_Q4_peaks$V2,
                   end = ED_Q4_peaks$V3),
)

mcols(gr_ED_Q1_peaks)$overlap <- countOverlaps(gr_ED_Q1_peaks, gr_Q1_ED_motifs)
mcols(gr_ED_Q2_peaks)$overlap <- countOverlaps(gr_ED_Q2_peaks, gr_Q2_ED_motifs)
mcols(gr_ED_Q3_peaks)$overlap <- countOverlaps(gr_ED_Q3_peaks, gr_Q3_ED_motifs)
mcols(gr_ED_Q4_peaks)$overlap <- countOverlaps(gr_ED_Q4_peaks, gr_Q4_ED_motifs)

#n of peaks with motifs
gr_ED_Q1_peaks[ gr_ED_Q1_peaks$overlap >= 1 ]
df_ED_Q1_peaks_overlapMotifs<- as.data.frame(gr_ED_Q1_peaks[ gr_ED_Q1_peaks$overlap >= 1 ])
nrow(df_ED_Q1_peaks_overlapMotifs)

gr_ED_Q2_peaks[ gr_ED_Q2_peaks$overlap >= 1 ]
df_ED_Q2_peaks_overlapMotifs<- as.data.frame(gr_ED_Q2_peaks[ gr_ED_Q2_peaks$overlap >= 1 ])
nrow(df_ED_Q2_peaks_overlapMotifs)

gr_ED_Q3_peaks[ gr_ED_Q3_peaks$overlap >= 1 ]
df_ED_Q3_peaks_overlapMotifs<- as.data.frame(gr_ED_Q3_peaks[ gr_ED_Q3_peaks$overlap >= 1 ])
nrow(df_ED_Q3_peaks_overlapMotifs)

gr_ED_Q4_peaks[ gr_ED_Q4_peaks$overlap >= 1 ]
df_ED_Q4_peaks_overlapMotifs<- as.data.frame(gr_ED_Q4_peaks[ gr_ED_Q4_peaks$overlap >= 1 ])
nrow(df_ED_Q4_peaks_overlapMotifs)


#prepare table for stats
ED_random <- cbind(df_random_ED_peaks_overlapMotifs, quartile= "1_random")
ED_Q1 <- cbind(df_ED_Q1_peaks_overlapMotifs, quartile= "2_Q1")
ED_Q2 <- cbind(df_ED_Q2_peaks_overlapMotifs, quartile= "3_Q2")
ED_Q3 <- cbind(df_ED_Q3_peaks_overlapMotifs, quartile= "4_Q3")
ED_Q4 <- cbind(df_ED_Q4_peaks_overlapMotifs, quartile= "4_Q4")
ED_all_stats <- rbind(ED_random, ED_Q1, ED_Q2, ED_Q3, ED_Q4)


# kruskal_test for multiple comparisons
res.kruskal <- ED_all_stats %>% rstatix::kruskal_test (overlap ~ quartile)
res.kruskal
#effect size
ED_all_stats %>% kruskal_effsize(overlap ~ quartile) #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# Pairwise comparisons using Dunn's test
pwc <- ED_all_stats %>% 
  dunn_test(overlap ~ quartile, p.adjust.method = "BH") 
pwc
fwrite(pwc, file ="MA02052/ED_DunnTest_BHpadj_nMotifs_per_peaks_counts_q0.05.tsv",sep= "\t", quote = F)

#WD violin plot score

ED_all_stats_violinPlot = ED_all_stats %>% group_by(quartile) %>% dplyr::summarise(num=n())

# Plot
Vplot_ED <- ED_all_stats %>%
  left_join(ED_all_stats_violinPlot) %>%
  mutate(myaxis = paste0(quartile, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=overlap, fill=quartile)) +
  geom_violin(width=1,alpha=0.5, color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 35)) +
  ylab("n Motifs overlap") + 
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("ED GAF peaks") +
  xlab("")
Vplot_ED
ggsave("MA02052/ED_GAF_quartiles_nMotifs.pdf", plot = Vplot_ED, width = 7, height = 5)

##################
#Same analysis for GA29 motif
### ED peaks motifs
ED_GAF_Q1_motifs <- read.table('GA29/ED_q1/fimo.tsv', header=T)
ED_GAF_Q2_motifs <- read.table('GA29/ED_q2/fimo.tsv', header=T)
ED_GAF_Q3_motifs <- read.table('GA29/ED_q3/fimo.tsv', header=T)
ED_GAF_Q4_motifs <- read.table('GA29/ED_q4/fimo.tsv', header=T)

#filter table by q-value <=0.05
ED_GAF_Q1_motifs <- ED_GAF_Q1_motifs[ED_GAF_Q1_motifs$q.value <= 0.05,]
ED_GAF_Q2_motifs <- ED_GAF_Q2_motifs[ED_GAF_Q2_motifs$q.value <= 0.05,]
ED_GAF_Q3_motifs <- ED_GAF_Q3_motifs[ED_GAF_Q3_motifs$q.value <= 0.05,]
ED_GAF_Q4_motifs <- ED_GAF_Q4_motifs[ED_GAF_Q4_motifs$q.value <= 0.05,]


ED_GAF_Q1_motifs <- ED_GAF_Q1_motifs[order(ED_GAF_Q1_motifs$sequence_name),]
ED_GAF_Q2_motifs <- ED_GAF_Q2_motifs[order(ED_GAF_Q2_motifs$sequence_name),]
ED_GAF_Q3_motifs <- ED_GAF_Q3_motifs[order(ED_GAF_Q3_motifs$sequence_name),]
ED_GAF_Q4_motifs <- ED_GAF_Q4_motifs[order(ED_GAF_Q4_motifs$sequence_name),]

# counting number of motifs
ED_GAF_Q1_motifs_count <- ED_GAF_Q1_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q2_motifs_count <- ED_GAF_Q2_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q3_motifs_count <- ED_GAF_Q3_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')
ED_GAF_Q4_motifs_count <- ED_GAF_Q4_motifs %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')


### random WD peaks
random_ED <- read.table('GA29/ED_random/fimo.tsv' , header=T)
random_ED <- random_ED[random_ED$q.value <= 0.05,]
random_ED <- random_ED[order(random_ED$sequence_name),]
random_ED_count <- random_ED %>% group_by(sequence_name) %>% summarise(total_count=n(),.groups = 'drop')


#stats score
#prepare table for stats
ED_random <- cbind(random_ED, class= "random")
ED_Q1 <- cbind(ED_GAF_Q1_motifs, class= "Q1")
ED_Q2 <- cbind(ED_GAF_Q2_motifs, class= "Q2")
ED_Q3 <- cbind(ED_GAF_Q3_motifs, class= "Q3")
ED_Q4 <- cbind(ED_GAF_Q4_motifs, class= "Q4")
ED_score_stats <- rbind(ED_random, ED_Q1, ED_Q2, ED_Q3, ED_Q4)


ED_score_stats %>% group_by(class) %>% get_summary_stats(score, type = "median_iqr")

# kruskal_test with Dunnets for multiple comparisons
ED_score_kruskal <- ED_score_stats %>% rstatix::kruskal_test (score ~ class)
ED_score_kruskal
#effect size
ED_score_stats %>% kruskal_effsize(score ~ class) #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# Pairwise comparisons using Dunn's test
pwc <- ED_score_stats %>% dunn_test(score ~ class, p.adjust.method = "BH") 
pwc
fwrite(pwc, file ="GA29/ED_DunnTest_BHpadj_scores_q005.tsv",sep= "\t", quote = F)
#end stats score

#calculate the total number of motifs in each class and create a vector.
ED_total_motifs <- c(sum(random_ED_count$total_count), sum(ED_GAF_Q1_motifs_count$total_count), sum(ED_GAF_Q2_motifs_count$total_count), sum(ED_GAF_Q3_motifs_count$total_count), sum(ED_GAF_Q4_motifs_count$total_count))
ED_total_motifs
names(ED_total_motifs) <- c("random", "GAF_Q1", "GAF_Q2", "GAF_Q3", "GAF_Q4")

# Convert to data frame
df <- data.frame(
  Group = factor(names(ED_total_motifs), levels = names(ED_total_motifs)),
  Count = as.numeric(ED_total_motifs)
)

v <- ggplot(df, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.9, direction = -1) +
  labs(title = "Eye disc", y = "Total GAF motifs", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(ylim = c(0, 35000))
v
ggsave("GA29/ED_Q1toQ4_loopAnchors_Total_n_motif_q0.05.pdf", plot = v, width = 7, height = 5)


#WD_random peaks with motifs
gr_random_ED_motifs <- GRanges(
  seqnames = random_ED$sequence_name,
  ranges = IRanges(start = random_ED$start,
                   end = random_ED$stop),
)

random_ED_peaks <- read.table('ED_random_bed_l967_n1610_seed24_NoChr.bed', header=F)
random_ED_peaks <- random_ED_peaks[order(random_ED_peaks$V1, random_ED_peaks$V2, random_ED_peaks$V3),]

gr_random_ED_peaks <- GRanges(
  seqnames = random_ED_peaks$V1,
  ranges = IRanges(start = random_ED_peaks$V2,
                   end = random_ED_peaks$V3),
)

overlap_random <- subsetByOverlaps(gr_random_ED_peaks, gr_random_ED_motifs)
mcols(gr_random_ED_peaks)$overlap <- countOverlaps(gr_random_ED_peaks, gr_random_ED_motifs)

#n of peaks with motifs
gr_random_ED_peaks[ gr_random_ED_peaks$overlap >= 1 ]
df_random_ED_peaks_overlapMotifs<- as.data.frame(gr_random_ED_peaks[ gr_random_ED_peaks$overlap >= 1 ])

# WD_Qn peaks with motifs
gr_Q1_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q1_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q1_motifs$start,
                   end = ED_GAF_Q1_motifs$stop),
)

gr_Q2_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q2_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q2_motifs$start,
                   end = ED_GAF_Q2_motifs$stop),
)

gr_Q3_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q3_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q3_motifs$start,
                   end = ED_GAF_Q3_motifs$stop),
)

gr_Q4_ED_motifs <- GRanges( 
  seqnames = ED_GAF_Q4_motifs$sequence_name,
  ranges = IRanges(start = ED_GAF_Q4_motifs$start,
                   end = ED_GAF_Q4_motifs$stop),
)

ED_Q1_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q1_NoChr.bed', header=F)
ED_Q1_peaks <- ED_Q1_peaks[order(ED_Q1_peaks$V1, ED_Q1_peaks$V2, ED_Q1_peaks$V3),]

ED_Q2_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q2_NoChr.bed', header=F)
ED_Q2_peaks <- ED_Q2_peaks[order(ED_Q2_peaks$V1, ED_Q2_peaks$V2, ED_Q2_peaks$V3),]

ED_Q3_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q3_NoChr.bed', header=F)
ED_Q3_peaks <- ED_Q3_peaks[order(ED_Q3_peaks$V1, ED_Q3_peaks$V2, ED_Q3_peaks$V3),]

ED_Q4_peaks <- read.table('ED_GAF_narrowPeaks_s100_FC2_q4_NoChr.bed', header=F)
ED_Q4_peaks <- ED_Q4_peaks[order(ED_Q4_peaks$V1, ED_Q4_peaks$V2, ED_Q4_peaks$V3),]

gr_ED_Q1_peaks <- GRanges( 
  seqnames = ED_Q1_peaks$V1,
  ranges = IRanges(start = ED_Q1_peaks$V2,
                   end = ED_Q1_peaks$V3),
)

gr_ED_Q2_peaks <- GRanges( 
  seqnames = ED_Q2_peaks$V1,
  ranges = IRanges(start = ED_Q2_peaks$V2,
                   end = ED_Q2_peaks$V3),
)

gr_ED_Q3_peaks <- GRanges( 
  seqnames = ED_Q3_peaks$V1,
  ranges = IRanges(start = ED_Q3_peaks$V2,
                   end = ED_Q3_peaks$V3),
)

gr_ED_Q4_peaks <- GRanges( 
  seqnames = ED_Q4_peaks$V1,
  ranges = IRanges(start = ED_Q4_peaks$V2,
                   end = ED_Q4_peaks$V3),
)

mcols(gr_ED_Q1_peaks)$overlap <- countOverlaps(gr_ED_Q1_peaks, gr_Q1_ED_motifs)
mcols(gr_ED_Q2_peaks)$overlap <- countOverlaps(gr_ED_Q2_peaks, gr_Q2_ED_motifs)
mcols(gr_ED_Q3_peaks)$overlap <- countOverlaps(gr_ED_Q3_peaks, gr_Q3_ED_motifs)
mcols(gr_ED_Q4_peaks)$overlap <- countOverlaps(gr_ED_Q4_peaks, gr_Q4_ED_motifs)

#n of peaks with motifs
gr_ED_Q1_peaks[ gr_ED_Q1_peaks$overlap >= 1 ]
df_ED_Q1_peaks_overlapMotifs<- as.data.frame(gr_ED_Q1_peaks[ gr_ED_Q1_peaks$overlap >= 1 ])
nrow(df_ED_Q1_peaks_overlapMotifs)

gr_ED_Q2_peaks[ gr_ED_Q2_peaks$overlap >= 1 ]
df_ED_Q2_peaks_overlapMotifs<- as.data.frame(gr_ED_Q2_peaks[ gr_ED_Q2_peaks$overlap >= 1 ])
nrow(df_ED_Q2_peaks_overlapMotifs)

gr_ED_Q3_peaks[ gr_ED_Q3_peaks$overlap >= 1 ]
df_ED_Q3_peaks_overlapMotifs<- as.data.frame(gr_ED_Q3_peaks[ gr_ED_Q3_peaks$overlap >= 1 ])
nrow(df_ED_Q3_peaks_overlapMotifs)

gr_ED_Q4_peaks[ gr_ED_Q4_peaks$overlap >= 1 ]
df_ED_Q4_peaks_overlapMotifs<- as.data.frame(gr_ED_Q4_peaks[ gr_ED_Q4_peaks$overlap >= 1 ])
nrow(df_ED_Q4_peaks_overlapMotifs)


#prepare table for stats
ED_random <- cbind(df_random_ED_peaks_overlapMotifs, quartile= "1_random")
ED_Q1 <- cbind(df_ED_Q1_peaks_overlapMotifs, quartile= "2_Q1")
ED_Q2 <- cbind(df_ED_Q2_peaks_overlapMotifs, quartile= "3_Q2")
ED_Q3 <- cbind(df_ED_Q3_peaks_overlapMotifs, quartile= "4_Q3")
ED_Q4 <- cbind(df_ED_Q4_peaks_overlapMotifs, quartile= "4_Q4")
ED_all_stats <- rbind(ED_random, ED_Q1, ED_Q2, ED_Q3, ED_Q4)


# kruskal_test for multiple comparisons
res.kruskal <- ED_all_stats %>% rstatix::kruskal_test (overlap ~ quartile)
res.kruskal
#effect size
ED_all_stats %>% kruskal_effsize(overlap ~ quartile) #The interpretation values commonly in published literature are: 0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
# Pairwise comparisons using Dunn's test
pwc <- ED_all_stats %>% 
  dunn_test(overlap ~ quartile, p.adjust.method = "BH") 
pwc
fwrite(pwc, file ="GA29/ED_DunnTest_BHpadj_nMotifs_per_peaks_counts_q0.05.tsv",sep= "\t", quote = F)

#WD violin plot score

ED_all_stats_violinPlot = ED_all_stats %>% group_by(quartile) %>% dplyr::summarise(num=n())

# Plot
Vplot_ED <- ED_all_stats %>%
  left_join(ED_all_stats_violinPlot) %>%
  mutate(myaxis = paste0(quartile, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=overlap, fill=quartile)) +
  geom_violin(width=1,alpha=0.5, color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 80)) +
  ylab("n Motifs overlap") + 
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("ED GAF peaks") +
  xlab("")
Vplot_ED
ggsave("GA29/ED_GAF_quartiles_nMotifs.pdf", plot = Vplot_ED, width = 7, height = 5)
