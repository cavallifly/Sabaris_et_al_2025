library(ggplot2)
library(viridis)
library(GenomicRanges)
library(dplyr)
library(tidyverse)
library(scales)

# setwd("D:/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/loop_distances_histogram/")
setwd("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/loop_distances_histogram/")

# Adjust the path and separator as needed
loops <- read.table("../mC_loops_Marco/Embryo.txt", header = FALSE, sep = "\t")
colnames(loops) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "score")

# Distance from end of anchor1 to start of anchor2
loops$distance <- loops$start2 - loops$start1
summary(loops)


# GAF peaks
gaf <- read.table("../bedtools/E_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
gaf_coords <- gaf[, c(1, 2, 3)]
colnames(gaf_coords) <- c("chr", "start", "end")

# GAF peaks
gr_gaf <- GRanges(seqnames = gaf_coords$chr,
                  ranges = IRanges(start = gaf_coords$start, end = gaf_coords$end))

# Loop anchors as GRanges
gr_anchor1 <- GRanges(seqnames = loops$chr1,
                      ranges = IRanges(start = loops$start1, end = loops$end1))
gr_anchor2 <- GRanges(seqnames = loops$chr2,
                      ranges = IRanges(start = loops$start2, end = loops$end2))

# Identify which loop anchors overlap GAF peaks
overlap1 <- subjectHits(findOverlaps(gr_anchor1, gr_gaf))
overlap2 <- subjectHits(findOverlaps(gr_anchor2, gr_gaf))

# Indexes of loops that overlap in at least one anchor
loop_hits <- unique(c(queryHits(findOverlaps(gr_anchor1, gr_gaf)),
                      queryHits(findOverlaps(gr_anchor2, gr_gaf))))

loops_filtered <- loops[loop_hits, ]
summary(loops_filtered)
write.table(loops_filtered, file = "NEW_E_filtered_loops_by_GAF_overlap.bedpe", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# proportion of loops =< 500 Kb
prop_gaf <- sum(loops_filtered$distance <= 500000) / nrow(loops_filtered)
prop_gaf

prop_all <- sum(loops$distance <= 500000) / nrow(loops)
prop_all

# Add a label column to distinguish the sets
loops$set <- "All loops"
loops_filtered$set <- "GAF-overlapping loops"

# Combine
combined_df <- rbind(
  loops[, c("distance", "set")],
  loops_filtered[, c("distance", "set")]
)


hist_E <- ggplot(combined_df, aes(x = distance, fill = set)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  coord_cartesian(xlim = c(0, 1000000)) +
  labs(
    title = "Loop size distribution (Embryo)",
    x = "Distance between loop anchors (bp)",
    y = "Density"
  ) +
  geom_vline(data = combined_df %>% group_by(set) %>% summarise(med = median(distance)),
             aes(xintercept = med, color = set),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("darkviolet", "yellow4")) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )
hist_E
ggsave("NEW_E_loopSizes_histogram_median.pdf", plot = hist_E, width = 7, height = 5)


################## Same than before for eye dics
# Adjust the path and separator as needed
ED_loops <- read.table("../mC_loops_Marco/ED.txt", header = T, sep = "\t")
colnames(ED_loops) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "score")

# Distance from end of anchor1 to start of anchor2
ED_loops$distance <- ED_loops$start2 - ED_loops$start1
summary(ED_loops)

# GAF intersected loop anchor file
ED_gaf <- read.table("../bedtools/ED_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
ED_gaf_coords <- ED_gaf[, c(1, 2, 3)]
colnames(ED_gaf_coords) <- c("chr", "start", "end")

# GAF peaks
ED_gr_gaf <- GRanges(seqnames = ED_gaf_coords$chr,
                     ranges = IRanges(start = ED_gaf_coords$start, end = ED_gaf_coords$end))

# Loop anchors as GRanges
ED_gr_anchor1 <- GRanges(seqnames = ED_loops$chr1,
                         ranges = IRanges(start = ED_loops$start1, end = ED_loops$end1))
ED_gr_anchor2 <- GRanges(seqnames = ED_loops$chr2,
                         ranges = IRanges(start = ED_loops$start2, end = ED_loops$end2))

# Identify which loop anchors overlap GAF peaks
ED_overlap1 <- subjectHits(findOverlaps(ED_gr_anchor1, ED_gr_gaf))
ED_overlap2 <- subjectHits(findOverlaps(ED_gr_anchor2, ED_gr_gaf))

# Indexes of loops that overlap in at least one anchor
ED_loop_hits <- unique(c(queryHits(findOverlaps(ED_gr_anchor1, ED_gr_gaf)),
                         queryHits(findOverlaps(ED_gr_anchor2, ED_gr_gaf))))

ED_loops_filtered <- ED_loops[ED_loop_hits, ]
summary(ED_loops_filtered)
write.table(ED_loops_filtered, file = "NEW_ED_filtered_loops_by_GAF_overlap.bedpe", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# proportion of loops =< 500 Kb
ED_prop_gaf <- sum(ED_loops_filtered$distance <= 500000) / nrow(ED_loops_filtered)
ED_prop_gaf

ED_prop_all <- sum(ED_loops$distance <= 500000) / nrow(ED_loops)
ED_prop_all

# Add a label column to distinguish the sets
ED_loops$set <- "All loops"
ED_loops_filtered$set <- "GAF-overlapping loops"

# Combine
ED_combined_df <- rbind(
  ED_loops[, c("distance", "set")],
  ED_loops_filtered[, c("distance", "set")]
)

#Histogram
hist_ED <-ggplot(ED_combined_df, aes(x = distance, fill = set)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  coord_cartesian(xlim = c(0, 1000000)) +
  labs(
    title = "Loop size distribution (Eye disc)",
    x = "Distance between loop anchors (bp)",
    y = "Density"
  ) +
  geom_vline(data = ED_combined_df %>% group_by(set) %>% summarise(med = median(distance)),
             aes(xintercept = med, color = set),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("darkviolet", "yellow4")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )
hist_ED
ggsave("NEW_ED_loopSizes_histogram_median.pdf", plot = hist_ED, width = 7, height = 5)

################## Same than before for wing dics
# Adjust the path and separator as needed
WD_loops <- read.table("../mC_loops_Marco/WD.txt", header = T, sep = "\t")
colnames(WD_loops) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "score")

# Distance from end of anchor1 to start of anchor2
WD_loops$distance <- WD_loops$start2 - WD_loops$start1
summary(WD_loops)

# GAF intersected loop anchor file
WD_gaf <- read.table("../bedtools/WD_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
WD_gaf_coords <- WD_gaf[, c(1, 2, 3)]
colnames(WD_gaf_coords) <- c("chr", "start", "end")

# GAF peaks
WD_gr_gaf <- GRanges(seqnames = WD_gaf_coords$chr,
                     ranges = IRanges(start = WD_gaf_coords$start, end = WD_gaf_coords$end))

# Loop anchors as GRanges
WD_gr_anchor1 <- GRanges(seqnames = WD_loops$chr1,
                         ranges = IRanges(start = WD_loops$start1, end = WD_loops$end1))
WD_gr_anchor2 <- GRanges(seqnames = WD_loops$chr2,
                         ranges = IRanges(start = WD_loops$start2, end = WD_loops$end2))

# Identify which loop anchors overlap GAF peaks
WD_overlap1 <- subjectHits(findOverlaps(WD_gr_anchor1, WD_gr_gaf))
WD_overlap2 <- subjectHits(findOverlaps(WD_gr_anchor2, WD_gr_gaf))

# Indexes of loops that overlap in at least one anchor
WD_loop_hits <- unique(c(queryHits(findOverlaps(WD_gr_anchor1, WD_gr_gaf)),
                         queryHits(findOverlaps(WD_gr_anchor2, WD_gr_gaf))))

WD_loops_filtered <- WD_loops[WD_loop_hits, ]
summary(WD_loops_filtered)
write.table(WD_loops_filtered, file = "NEW_WD_filtered_loops_by_GAF_overlap.bedpe", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# proportion of loops =< 500 Kb
WD_prop_gaf <- sum(WD_loops_filtered$distance <= 500000) / nrow(WD_loops_filtered)
WD_prop_gaf

WD_prop_all <- sum(WD_loops$distance <= 500000) / nrow(WD_loops)
WD_prop_all

# Add a label column to distinguish the sets
WD_loops$set <- "All loops"
WD_loops_filtered$set <- "GAF-overlapping loops"

# Combine
WD_combined_df <- rbind(
  WD_loops[, c("distance", "set")],
  WD_loops_filtered[, c("distance", "set")]
)

#Histogram
hist_WD <-ggplot(WD_combined_df, aes(x = distance, fill = set)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  coord_cartesian(xlim = c(0, 1000000)) +
  labs(
    title = "Loop size distribution (Wing disc)",
    x = "Distance between loop anchors (bp)",
    y = "Density"
  ) +
  geom_vline(data = WD_combined_df %>% group_by(set) %>% summarise(med = median(distance)),
             aes(xintercept = med, color = set),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("darkviolet", "yellow4")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )
hist_WD
ggsave("NEW_WD_loopSizes_histogram_median.pdf", plot = hist_WD, width = 7, height = 5)


########## barPlot for Loops overlapping GAF peaks (at any anchor)
E_total_anchors <- nrow(loops)
E_n_overlap <- nrow(loops_filtered)
E_n_nonoverlap <- E_total_anchors - E_n_overlap

E_df_summary <- data.frame(
  Condition = "Embryo",
  Category = c("Overlap", "No Overlap"),
  Count = c(E_n_overlap, E_n_nonoverlap)
)
E_df_summary$Proportion <- E_df_summary$Count / sum(E_df_summary$Count)


ED_total_anchors <- nrow(ED_loops)
ED_n_overlap <- nrow(ED_loops_filtered)
ED_n_nonoverlap <- ED_total_anchors - ED_n_overlap

ED_df_summary <- data.frame(
  Condition = "Eye disc",
  Category = c("Overlap", "No Overlap"),
  Count = c(ED_n_overlap, ED_n_nonoverlap)
)
ED_df_summary$Proportion <- ED_df_summary$Count / sum(ED_df_summary$Count)

WD_total_anchors <- nrow(WD_loops)
WD_n_overlap <- nrow(WD_loops_filtered)
WD_n_nonoverlap <- WD_total_anchors - WD_n_overlap

WD_df_summary <- data.frame(
  Condition = "Wing disc",
  Category = c("Overlap", "No Overlap"),
  Count = c(WD_n_overlap, WD_n_nonoverlap)
)
WD_df_summary$Proportion <- WD_df_summary$Count / sum(WD_df_summary$Count)

#---------------------------
# Combine and plot
#---------------------------
df_summary <- rbind(E_df_summary, ED_df_summary, WD_df_summary)

pdf("NEW_stackedBarPlot_E_ED_WD_Loops_intersected_GAFpeaks.pdf")
ggplot(df_summary, aes(x = Condition, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Overlap" = "#414487FF", "No Overlap" = "#22A884FF")) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(title = "Loops overlapping GAF peaks",
       x = "", y = "Proportion of Peaks") +
  theme_minimal(base_size = 18) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()

