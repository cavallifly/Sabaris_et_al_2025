### What do you need to run the analysis and obtain the figures? (AKA Dependencies) ###
#misha
#shaman
#ggpubr
#dplyr
#doParallel
#bash: identify, convert

scriptsDir=${PWD}/scripts_clean/

# list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed
samples="larvae_DWT larvae_DPRE1 larvae_double LD_gypsy2"

# Generate the modified dm6 reference genome
#bash     ${scriptsDir}/01_editFasta.sh &> 01_editFasta.out

# Mapping of the datasets and the mishaDB entries for the observed tracks
# have been obtained using the "shHiC2" pipeline in https://github.com/tanaylab/schic2
# T. Nagano, Y. Lubling, C. Várnai, C. Dudley, W. Leung, Y. Baran, N. M. Cohen, S. Wingett, P. Fraser, A. Tanay,
# Cell-cycle dynamics of chromosomal organization at single-cell resolution. Nature 547, 61–67 (2017).

# Compute statistics on the mapped samples
#bash     ${scriptsDir}/01_computeValidPairsFromTracks.sh &> 01_computeValidPairsFromTracks.out
#bash     ${scriptsDir}/01_count_readPairs_from_fastq.sh  &> 01_count_readPairs_from_fastq.out

# Compute the shuffled (expected) tracks and the shaman Hi-C scores
# The misha tracks for the observed counts must be downloaded from following
# GEO entries:
#GSM7888066	Larv WT HiC Repli1
#GSM7888067	Larv WT HiC Repli2
#GSM7888068	Larv WT HiC Repli3
#GSM7888077	Larv gypsy2 insertion HiC Repli1
#GSM7888078	Larv gypsy2 insertion HiC Repli2
#GSM7888079	Larv gypsy2 insertion HiC Repli3
# ADD missing GEO accession codes XXX

#Rscript ${scriptsDir}/02_generateShuffleTrack.R &>> 02_generateShuffleTrack.out
#Rscript ${scriptsDir}/03_computeScoreTrack.R    &>> 03_computeScoreTrack.out

# Analysis of the Hi-C maps
#Rscript ${scriptsDir}/04_computeInsulationTracks.R        &>> 04_computeInsulationTracks.out
#Rscript ${scriptsDir}/05_computeDownSampledTracks.R       &>> 05_computeDownSampledTracks.out
#Rscript ${scriptsDir}/06_computeDiffScoreTrack_parallel.R &>> 06_computeDiffScoreTrack_parallel.out

# Data and Figures for each figure panel
Rscript ${scriptsDir}/07_ScoreMapsPlot_Fig1C.R                 &>> 07_ScoreMapsPlot_Fig1C.out
Rscript ${scriptsDir}/08_INSplot_Fig1D.R                       &>> 08_INSplot_Fig1D.out
Rscript ${scriptsDir}/09_INSquantification_Fig1E.R             &>> 09_INSquantification_Fig1E.out
Rscript ${scriptsDir}/10_DiffScoreMapsPlot_Fig1F.R             &>> 10_DiffScoreMapsPlot_Fig1F.out
Rscript ${scriptsDir}/11_ScoreMapsQuantification_Fig1G.R       &>> 11_ScoreMapsQuantification_Fig1G.out
Rscript ${scriptsDir}/12a_get_obsContacts_PcG_domains_Fig1H.R  &>> 12a_get_obsContacts_PcG_domains_Fig1H.out
conda run -n DEseq2 Rscript ${scriptsDir}/12b_plot_obsContacts_PcG_domains_Fig1H.R &>> 12b_plot_obsContacts_PcG_domains_Fig1H.out
Rscript ${scriptsDir}/13_ScoreMapsPlot_Fig2C.R                 &>> 13_ScoreMapsPlot_Fig2C.out
Rscript ${scriptsDir}/14_INSplot_Fig2D.R                       &>> 14_INSplot_Fig2D.out
Rscript ${scriptsDir}/15_INSquantification_Fig2E.R             &>> 15_INSquantification_Fig2E.out
Rscript ${scriptsDir}/16_DiffScoreMapsPlot_Fig2G.R             &>> 16_DiffScoreMapsPlot_Fig2G.out
Rscript ${scriptsDir}/17a_ScoreMapsQuantification_Fig2Fa.R     &>>  17a_ScoreMapsQuantification_Fig2Fa.out
Rscript ${scriptsDir}/17b_ScoreMapsQuantification_Fig2Fb.R     &>>  17b_ScoreMapsQuantification_Fig2Fb.out

# Obtain the triangular maps
bash ${scriptsDir}/18_generate_triangular_pngs.sh &> 18_generate_triangular_pngs.out

# Move results in folders
mkdir -p Data_for_figures
mkdir -p Figure_panels

mv *.png *.pdf Figure_panels
mv *.tab *.tsv Data_for_figures/
