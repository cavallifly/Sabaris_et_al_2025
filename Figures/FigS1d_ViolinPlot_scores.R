# Load required packages
library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(rstatix)
library(ggpubr)
library(data.table)


setwd("D:/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/")

## Embryo peaks
# Load GAF peaks and loop overlaps (adjust paths as needed)
E_GAF_npeaks <- read.delim("E_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
E_GAF_npeaks_inters_E_loops <- read.delim("E_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")

# Deduplicate intersected peaks
vector_E_GAF_loops <- unique(as.character(E_GAF_npeaks_inters_E_loops$V4))
E_GAF_npeaks$V4 <- as.character(E_GAF_npeaks$V4)

# Label classes
E_all <- E_GAF_npeaks %>%
  mutate(class = "all")

E_loop <- E_GAF_npeaks %>%
  filter(V4 %in% vector_E_GAF_loops) %>%
  mutate(class = "loop")

E_NO_loop <- E_GAF_npeaks %>%
  filter(!V4 %in% vector_E_GAF_loops) %>%
  mutate(class = "NO_loop")

# Combine all into one dataset
E_score_stats <- bind_rows(E_all, E_NO_loop, E_loop)

# Rename score column (assumes V5 holds the peak score)
colnames(E_score_stats)[5] <- "Score"

# Summary stats (optional)
E_score_stats %>% group_by(class) %>% get_summary_stats(Score, type = "median_iqr")
E_score_stats %>% kruskal_test(Score ~ class)
pwc <- E_score_stats %>% dunn_test(Score ~ class, p.adjust.method = "BH")
pwc
fwrite(pwc, file ="E_dunnTest_BHpadj_scores_loops.tsv",sep= "\t", quote = F)

# Count per class for n labels
E_score_value <- E_score_stats %>% 
  group_by(class) %>% 
  summarise(num = n())

# Join and create x-axis label with sample size
E_score_plot <- E_score_stats %>%
  left_join(E_score_value, by = "class") %>%
  mutate(myaxis = paste0(class, "\n", "n=", num))

# Plot
Vplot <- ggplot(E_score_plot, aes(x = myaxis, y = Score, fill = class)) +
  geom_violin(width = 1, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 6000)) +
  ylab("Peak score") +
  xlab("") +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Embryo GAF peaks")
Vplot

# Save plot
ggsave("E_violinPlot_E_GAF_narrowPeaks_score_loops.pdf", plot = Vplot, width = 7, height = 5)

## WD peaks
# Load GAF peaks and loop overlaps (adjust paths as needed)
WD_GAF_npeaks <- read.delim("WD_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
WD_GAF_npeaks_inters_E_loops <- read.delim("WD_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")

# Deduplicate intersected peaks
vector_WD_GAF_loops <- unique(as.character(WD_GAF_npeaks_inters_E_loops$V4))
WD_GAF_npeaks$V4 <- as.character(WD_GAF_npeaks$V4)

# Label classes
WD_all <- WD_GAF_npeaks %>%
  mutate(class = "all")

WD_loop <- WD_GAF_npeaks %>%
  filter(V4 %in% vector_WD_GAF_loops) %>%
  mutate(class = "loop")

WD_NO_loop <- WD_GAF_npeaks %>%
  filter(!V4 %in% vector_WD_GAF_loops) %>%
  mutate(class = "NO_loop")

# Combine all into one dataset
WD_score_stats <- bind_rows(WD_all, WD_NO_loop, WD_loop)

# Rename score column (assumes V5 holds the peak score)
colnames(WD_score_stats)[5] <- "Score"

# Summary stats (optional)
WD_score_stats %>% group_by(class) %>% get_summary_stats(Score, type = "median_iqr")
WD_score_stats %>% kruskal_test(Score ~ class)
WD_pwc <- WD_score_stats %>% dunn_test(Score ~ class, p.adjust.method = "BH")
WD_pwc
fwrite(WD_pwc, file ="WD_dunnTest_BHpadj_scores_loops.tsv",sep= "\t", quote = F)

# Count per class for n labels
WD_score_value <- WD_score_stats %>% 
  group_by(class) %>% 
  summarise(num = n())

# Join and create x-axis label with sample size
WD_score_plot <- WD_score_stats %>%
  left_join(WD_score_value, by = "class") %>%
  mutate(myaxis = paste0(class, "\n", "n=", num))

# Plot
WD_Vplot <- ggplot(WD_score_plot, aes(x = myaxis, y = Score, fill = class)) +
  geom_violin(width = 1, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 10000)) +
  ylab("Peak score") +
  xlab("") +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Wing disc GAF peaks")
WD_Vplot

# Save plot
ggsave("WD_violinPlot_E_GAF_narrowPeaks_score_loops.pdf", plot = WD_Vplot, width = 7, height = 5)


## ED peaks
# Load GAF peaks and loop overlaps (adjust paths as needed)
ED_GAF_npeaks <- read.delim("ED_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep = "\t")
ED_GAF_npeaks_inters_ED_loops <- read.delim("ED_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep = "\t")

# Deduplicate intersected peaks
vector_ED_GAF_loops <- unique(as.character(ED_GAF_npeaks_inters_ED_loops$V4))
ED_GAF_npeaks$V4 <- as.character(ED_GAF_npeaks$V4)

# Label classes
ED_all <- ED_GAF_npeaks %>%
  mutate(class = "all")

ED_loop <- ED_GAF_npeaks %>%
  filter(V4 %in% vector_ED_GAF_loops) %>%
  mutate(class = "loop")

ED_NO_loop <- ED_GAF_npeaks %>%
  filter(!V4 %in% vector_ED_GAF_loops) %>%
  mutate(class = "NO_loop")

# Combine all into one dataset
ED_score_stats <- bind_rows(ED_all, ED_NO_loop, ED_loop)

# Rename score column (assumes V5 holds the peak score)
colnames(ED_score_stats)[5] <- "Score"

# Summary stats (optional)
ED_score_stats %>% group_by(class) %>% get_summary_stats(Score, type = "median_iqr")
ED_score_stats %>% kruskal_test(Score ~ class)
ED_pwc <- ED_score_stats %>% dunn_test(Score ~ class, p.adjust.method = "BH")
ED_pwc
fwrite(ED_pwc, file ="ED_dunnTest_BHpadj_scores_loops.tsv",sep= "\t", quote = F)

# Count per class for n labels
ED_score_value <- ED_score_stats %>% 
  group_by(class) %>% 
  summarise(num = n())

# Join and create x-axis label with sample size
ED_score_plot <- ED_score_stats %>%
  left_join(ED_score_value, by = "class") %>%
  mutate(myaxis = paste0(class, "\n", "n=", num))

# Plot
ED_Vplot <- ggplot(ED_score_plot, aes(x = myaxis, y = Score, fill = class)) +
  geom_violin(width = 1, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 10000)) +
  ylab("Peak score") +
  xlab("") +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Eye disc GAF peaks")
ED_Vplot

# Save plot
ggsave("ED_violinPlot_ED_GAF_narrowPeaks_score_loops.pdf", plot = ED_Vplot, width = 7, height = 5)
