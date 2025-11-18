##########################################################
## LOOP–GAF OVERLAP PIPELINE
## - Anchor expansion
## - Loop classification
## - Per-class BEDPE export
## - Summary tables
## - Stacked bar plot
##########################################################

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)

##########################################################
# FUNCTION: Run full analysis for one loop + GAF set
##########################################################

analyze_loops_gaf <- function(loop_file,
                              gaf_file,
                              condition_name,
                              # anchor_extension = 2000,
                              output_prefix = condition_name) {
  
  message("Processing condition: ", condition_name)
  
  # ------------------------------------------------------------
  # 1. Load loops (BEDPE)
  # ------------------------------------------------------------
  loops <- read.table(loop_file, header = FALSE)
  colnames(loops)[1:6] <- c("chr1","start1","end1","chr2","start2","end2")
  
  # # Expand anchors ±2000 bp (for 4kb windows)
  # loops$start1 <- pmax(0, loops$start1 - anchor_extension)
  # loops$end1   <- loops$end1 + anchor_extension
  # loops$start2 <- pmax(0, loops$start2 - anchor_extension)
  # loops$end2   <- loops$end2 + anchor_extension
  
  # Convert to GRanges
  gr_anchor1 <- GRanges(loops$chr1, IRanges(loops$start1, loops$end1))
  gr_anchor2 <- GRanges(loops$chr2, IRanges(loops$start2, loops$end2))
  
  # ------------------------------------------------------------
  # 2. Read full GAF peak file
  # ------------------------------------------------------------
  gaf <- read.table(gaf_file, header = FALSE)
  gaf <- gaf[,1:3]
  colnames(gaf) <- c("chr","start","end")
  
  gr_gaf <- GRanges(gaf$chr, IRanges(gaf$start, gaf$end))
  
  # ------------------------------------------------------------
  # 3. Compute overlaps
  # ------------------------------------------------------------
  ov1 <- queryHits(findOverlaps(gr_anchor1, gr_gaf))
  ov2 <- queryHits(findOverlaps(gr_anchor2, gr_gaf))
  
  # Loop-level classification
  nloops <- nrow(loops)
  classification <- rep("none", nloops)
  
  classification[unique(ov1)] <- "anchor1_only"
  classification[unique(ov2)] <- ifelse(
    classification[unique(ov2)] == "anchor1_only",
    "both",
    "anchor2_only"
  )
  
  # Collapsed classification
  collapsed <- dplyr::recode(classification,
                             anchor1_only = "one_anchor",
                             anchor2_only = "one_anchor")
  
  # ------------------------------------------------------------
  # 4. Add distance between loop anchors
  # ------------------------------------------------------------
  loops$distance <- loops$start2 - loops$start1
  
  # ------------------------------------------------------------
  # 5. Build final output table
  # ------------------------------------------------------------
  df <- loops %>%
    mutate(classification = classification,
           collapsed_class = collapsed,
           condition = condition_name)
  
  out_tsv <- paste0(output_prefix, "_loops_with_GAF_classification.tsv")
  write.table(df, out_tsv, sep="\t", quote=FALSE, row.names=FALSE)
  message(" → Saved: ", out_tsv)
  
  # ------------------------------------------------------------
  # 6. Export per-category BEDPE
  # ------------------------------------------------------------
  save_cat <- function(cat) {
    fn <- paste0(output_prefix, "_", cat, ".bedpe")
    write.table(df[df$collapsed_class == cat, 1:6],
                fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    message(" → Saved: ", fn)
  }
  
  save_cat("none")
  save_cat("one_anchor")
  save_cat("both")
  
  # ------------------------------------------------------------
  # 7. Summary table
  # ------------------------------------------------------------
  summary_df <- df %>%
    group_by(collapsed_class) %>%
    summarise(n = n(),
              proportion = n / nloops) %>%
    mutate(condition = condition_name)
  
  return(summary_df)
}

##########################################################
# FUNCTION: Create stacked bar plot for many conditions
##########################################################
setwd("/Users/gonzalosabaris/Nextcloud/Gonzalez_lab/Writting_papers/2024_GAF_paper/2025/2025_07_Analysis_paper/Overlap_GAF_loops/")

plot_loop_overlap <- function(summary_df, outfile="loop_GAF_overlap_stackedBars.pdf") {
  
  summary_df$collapsed_class <- factor(summary_df$collapsed_class,
                                       levels=c("none","one_anchor","both"))
  
  p <- ggplot(summary_df,
              aes(x=condition, y=proportion, fill=collapsed_class)) +
    geom_bar(stat="identity", width=0.6, color="black") +
    geom_text(aes(label=scales::percent(proportion, accuracy=0.1)),
              position=position_stack(vjust=0.5),
              color="white", size=4) +
    scale_fill_manual(values = c(
      none        = "#22A884FF",  # green
      one_anchor  = "#414487FF",  # blue
      both        = "#FDE725FF"   # yellow
    )) +
    scale_y_continuous(labels=percent_format()) +
    labs(title="Loops overlapping GAF peaks",
         x="", y="Proportion of loops", fill="Category") +
    theme_minimal(base_size=16) +
    theme(plot.title = element_text(hjust=0.5))
  
  ggsave(outfile, p, width=8, height=6)
  message(" → Saved plot: ", outfile)
}

##########################################################
# MULTI-DATASET DRIVER (EDIT THIS PART)
##########################################################

# Example:
# Embryo, Eye disc, Wing disc
loop_files <- c(
  Embryo="../mC_loops_Marco/Embryo.txt",
  Eye="../mC_loops_Marco/ED.txt",
  Wing="../mC_loops_Marco/WD.txt"
)

gaf_files <- c(
  Embryo="../bedtools/E_GAF_narrowPeaks_s100_FC2.narrowPeak",
  Eye="../bedtools/ED_GAF_narrowPeaks_s100_FC2.narrowPeak",
  Wing="../bedtools/WD_GAF_narrowPeaks_s100_FC2.narrowPeak"
)

# Run all
all_summary <- data.frame()

for(cond in names(loop_files)) {
  s <- analyze_loops_gaf(loop_files[[cond]],
                         gaf_files[[cond]],
                         condition_name = cond,
                         # anchor_extension = 2000,
                         output_prefix = paste0(cond, "_GAF"))
  all_summary <- rbind(all_summary, s)
}

# Plot combined results
plot_loop_overlap(all_summary,
                  outfile = "stackedBarPlot_allConditions_loops_GAF.pdf")

##########################################################
message("DONE ✔✔✔")
