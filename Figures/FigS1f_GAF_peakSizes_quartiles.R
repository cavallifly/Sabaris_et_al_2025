require(data.table)
library(ggplot2)
library(hrbrthemes)
setwd("D:/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/")

#peak length distribution per quartiles
E_GAF_npeaks <- read.delim("E_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
E_GAF_npeaks$peak_size <- E_GAF_npeaks$V3 - E_GAF_npeaks$V2
head(E_GAF_npeaks)

E_GAF_npeaks$V11 <- as.factor(E_GAF_npeaks$V11)
levels(E_GAF_npeaks$V11) <- c("Q1", "Q2", "Q3", "Q4")  # Rename quartiles

# Basic violin plot
p <- ggplot(E_GAF_npeaks, aes(x=V11, y=peak_size, fill=V11)) + 
  geom_violin(width=1,alpha=0.5, color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 6000)) +
  ylab("Peak size (bp)") + 
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Embryo GAF peaks (quartiles)") +
  xlab("")
p
ggsave("E_GAF_quartiles_bySize.pdf", plot = p, width = 7, height = 5)



WD_GAF_npeaks <- read.delim("WD_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
WD_GAF_npeaks$peak_size <- WD_GAF_npeaks$V3 - WD_GAF_npeaks$V2
head(WD_GAF_npeaks)

WD_GAF_npeaks$V11 <- as.factor(WD_GAF_npeaks$V11)
levels(WD_GAF_npeaks$V11) <- c("Q1", "Q2", "Q3", "Q4")  # Rename quartiles

p_WD <- ggplot(WD_GAF_npeaks, aes(x=V11, y=peak_size, fill=V11)) + 
  geom_violin(width=1,alpha=0.5, color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 6000)) +
  ylab("Peak size (bp)") + 
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Wing disc GAF peaks (quartiles)") +
  xlab("")
p_WD
ggsave("WD_GAF_quartiles_bySize.pdf", plot = p_WD, width = 7, height = 5)


## ED
ED_GAF_npeaks <- read.delim("ED_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
ED_GAF_npeaks$peak_size <- ED_GAF_npeaks$V3 - ED_GAF_npeaks$V2
head(ED_GAF_npeaks)

ED_GAF_npeaks$V11 <- as.factor(ED_GAF_npeaks$V11)
levels(ED_GAF_npeaks$V11) <- c("Q1", "Q2", "Q3", "Q4")  # Rename quartiles

p_ED <- ggplot(ED_GAF_npeaks, aes(x=V11, y=peak_size, fill=V11)) + 
  geom_violin(width=1,alpha=0.5, color = NA) +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 6000)) +
  ylab("Peak size (bp)") + 
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_ipsum(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Eye disc GAF peaks (quartiles)") +
  xlab("")
p_ED
ggsave("ED_GAF_quartiles_bySize.pdf", plot = p_ED, width = 7, height = 5)
