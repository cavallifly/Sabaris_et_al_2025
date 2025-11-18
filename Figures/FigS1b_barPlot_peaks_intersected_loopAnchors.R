library(ggplot2)
library(scales)

setwd("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/")
# setwd("/Users/gonzalosabaris/Nextcloud/NGS_Data/synPRE_CnR/CnR_GAF_triplicate_LF/macs3")
# setwd("D:/Nextcloud/NGS_Data/synPRE_CnR/CnR_GAF_triplicate_LF/macs3")

##load peaks
#WD
WD_GAF_npeaks <- read.delim("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/WD_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
WD_GAF_npeaks_inters_WD_loops <- read.delim("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/WD_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep="\t")

#Embryo
E_GAF_npeaks <- read.delim("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/E_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
E_GAF_npeaks_inters_E_loops <- read.delim("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/bedtools/E_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep="\t")

#ED
ED_GAF_npeaks <- read.delim("ED_GAF_narrowPeaks_s100_FC2.narrowPeak", header = FALSE, sep="\t")
ED_GAF_npeaks_inters_ED_loops <- read.delim("ED_GAFpeaks_intersected_anchLoops.bed", header = FALSE, sep="\t")

#---------------------------
# Load and process Embryo
#---------------------------
E_all_peaks <- unique(as.character(E_GAF_npeaks$V4))
E_overlapping_peaks <- unique(as.character(E_GAF_npeaks_inters_E_loops$V4))

E_total_peaks <- length(E_all_peaks)
E_n_overlap <- length(E_overlapping_peaks)
E_n_nonoverlap <- E_total_peaks - E_n_overlap

E_df_summary <- data.frame(
  Condition = "Embryo",
  Category = c("Overlap", "No Overlap"),
  Count = c(E_n_overlap, E_n_nonoverlap)
)
E_df_summary$Proportion <- E_df_summary$Count / sum(E_df_summary$Count)

#---------------------------
# Load and process WD
#---------------------------
WD_all_peaks <- unique(as.character(WD_GAF_npeaks$V4))
WD_overlapping_peaks <- unique(as.character(WD_GAF_npeaks_inters_WD_loops$V4))

WD_total_peaks <- length(WD_all_peaks)
WD_n_overlap <- length(WD_overlapping_peaks)
WD_n_nonoverlap <- WD_total_peaks - WD_n_overlap

WD_df_summary <- data.frame(
  Condition = "Wing disc",
  Category = c("Overlap", "No Overlap"),
  Count = c(WD_n_overlap, WD_n_nonoverlap)
)
WD_df_summary$Proportion <- WD_df_summary$Count / sum(WD_df_summary$Count)

#---------------------------
# Load and process ED
#---------------------------
ED_all_peaks <- unique(as.character(ED_GAF_npeaks$V4))
ED_overlapping_peaks <- unique(as.character(ED_GAF_npeaks_inters_ED_loops$V4))

ED_total_peaks <- length(ED_all_peaks)
ED_n_overlap <- length(ED_overlapping_peaks)
ED_n_nonoverlap <- ED_total_peaks - ED_n_overlap

ED_df_summary <- data.frame(
  Condition = "Eye disc",
  Category = c("Overlap", "No Overlap"),
  Count = c(ED_n_overlap, ED_n_nonoverlap)
)
ED_df_summary$Proportion <- ED_df_summary$Count / sum(ED_df_summary$Count)

#---------------------------
# Combine and plot
#---------------------------
df_summary <- rbind(E_df_summary, WD_df_summary, ED_df_summary)

pdf("stackedBarPlot_E_WD_ED_peaks_intersected_loopAnchors.pdf")
ggplot(df_summary, aes(x = Condition, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Overlap" = "#414487FF", "No Overlap" = "#22A884FF")) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  geom_text(aes(label = percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),  # put labels inside bars
            color = "white", size = 5) +
  labs(title = "GAF peaks overlapping loop anchors",
       x = "", y = "Proportion of Peaks") +
  theme_minimal(base_size = 18) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

dev.off()
