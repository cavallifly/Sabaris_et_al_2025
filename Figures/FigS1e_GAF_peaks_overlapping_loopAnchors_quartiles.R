# loopAnchors intersection with quartiles.
library(dplyr)
library(ggplot2)
library(viridis)
setwd("D:/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/")

# -------------------------
# Embryo: count overlaps per quartile
# -------------------------
E_GAF_npeaks_inters_E_loops <- read.delim("E_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")
E_GAF_npeaks_inters_E_loops <- subset(E_GAF_npeaks_inters_E_loops, select = c(1:11))
E_GAF_npeaks_inters_E_loops_unique <- distinct(E_GAF_npeaks_inters_E_loops, V4, .keep_all = TRUE)

E_GAF_total_q1 <- sum(E_GAF_npeaks_inters_E_loops_unique$V11 == '1')
E_GAF_total_q2 <- sum(E_GAF_npeaks_inters_E_loops_unique$V11 == '2')
E_GAF_total_q3 <- sum(E_GAF_npeaks_inters_E_loops_unique$V11 == '3')
E_GAF_total_q4 <- sum(E_GAF_npeaks_inters_E_loops_unique$V11 == '4')

# -------------------------
# WD: count overlaps per quartile
# -------------------------
WD_GAF_npeaks_inters_WD_loops <- read.delim("WD_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")
WD_GAF_npeaks_inters_WD_loops <- subset(WD_GAF_npeaks_inters_WD_loops, select = c(1:11))
WD_GAF_npeaks_inters_WD_loops_unique <- distinct(WD_GAF_npeaks_inters_WD_loops, V4, .keep_all = TRUE)

WD_GAF_total_q1 <- sum(WD_GAF_npeaks_inters_WD_loops_unique$V11 == '1')
WD_GAF_total_q2 <- sum(WD_GAF_npeaks_inters_WD_loops_unique$V11 == '2')
WD_GAF_total_q3 <- sum(WD_GAF_npeaks_inters_WD_loops_unique$V11 == '3')
WD_GAF_total_q4 <- sum(WD_GAF_npeaks_inters_WD_loops_unique$V11 == '4')

# -------------------------
# ED: count overlaps per quartile
# -------------------------
ED_GAF_npeaks_inters_ED_loops <- read.delim("ED_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")
ED_GAF_npeaks_inters_ED_loops <- subset(ED_GAF_npeaks_inters_ED_loops, select = c(1:11))
ED_GAF_npeaks_inters_ED_loops_unique <- distinct(ED_GAF_npeaks_inters_ED_loops, V4, .keep_all = TRUE)

ED_GAF_total_q1 <- sum(ED_GAF_npeaks_inters_ED_loops_unique$V11 == '1')
ED_GAF_total_q2 <- sum(ED_GAF_npeaks_inters_ED_loops_unique$V11 == '2')
ED_GAF_total_q3 <- sum(ED_GAF_npeaks_inters_ED_loops_unique$V11 == '3')
ED_GAF_total_q4 <- sum(ED_GAF_npeaks_inters_ED_loops_unique$V11 == '4')


# -------------------------
# Create summary data frame
# -------------------------
df <- data.frame(
  Condition = rep(c("Embryo", "ED", "WD"), each = 4),
  Quartile = rep(c("Q4", "Q3", "Q2", "Q1"), 3),
  OverlapCount = c(E_GAF_total_q4, E_GAF_total_q3, E_GAF_total_q2, E_GAF_total_q1,
                   ED_GAF_total_q4, ED_GAF_total_q3, ED_GAF_total_q2, ED_GAF_total_q1,
                   WD_GAF_total_q4, WD_GAF_total_q3, WD_GAF_total_q2, WD_GAF_total_q1)
)

# Normalize to proportions within each condition
df <- df %>%
  group_by(Condition) %>%
  mutate(Proportion = OverlapCount / sum(OverlapCount))

# -------------------------
# Plot stacked bar plot
# -------------------------
barPlot <- ggplot(df, aes(x = Condition, y = Proportion, fill = Quartile)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.9, direction = -1) +  
  labs(title = "GAF peaks overlapping loop anchors by quartile",
       x = "", y = "Proportion of Peaks") +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
barPlot
# Save plot
ggsave("NEW_GAF_peaks_overlapping_loopAnchors_quartiles.pdf", plot = barPlot, width = 7, height = 5)
