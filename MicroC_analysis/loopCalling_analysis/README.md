# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

Ensure that you have a working version of R. The scripts in this repository have been tested for version 4.0.5 (2021-03-31).
<!--
Install the following R packages: ggplot2 and DiffBind.
-->
The programs pyGenomeTracks, samtools, sambamba, deeptools, and bowtie2.

# Input data #
Next, you should download the .mcool files from the GEO entry GSE310299 and place them in the mcoolFiles folder.

Now, you are ready to run the scripts in ./scripts.

and run one after the other the following commands following the messages:
```
scriptsDir=./scripts/
bash ${scriptsDir}/01_classify_Furlong_loops.sh
bash ${scriptsDir}/01_plot_looping_regions_Furlong.cmd
bash ${scriptsDir}/02_callLoops_using_mustache_launch.sh
bash ${scriptsDir}/03_remove_periCentromeric_regions.sh
bash ${scriptsDir}/04_filter_loops_by_FDRperResolution.sh
bash ${scriptsDir}/05_check_Nloops_per_FDR_threshold_files.sh
bash ${scriptsDir}/05_rank_loopDetection_by_FDRperResolution.sh
bash ${scriptsDir}/05_plot_looping_regions_for_visualScreening.sh
bash ${scriptsDir}/06_get_accepted_loops_loops.sh
bash ${scriptsDir}/07_recenter_loops_by_maximum_shamanScore.R
bash ${scriptsDir}/07_recenter_loops_by_maximum_shamanScore.sh
```

Pay attention that some part involve visual screening of the loops in the generated image files.

Once the procedure is finished, you will obtain the list of loops per condition.

-->

## Contributions ##
The code in this repository has been developed at the [Cavalli Lab](https://www.igh.cnrs.fr/en/research/departments/genome-dynamics/chromatin-and-cell-biology) with the contributions of Marco Di Stefano, Gonzalo Sabaris, and Giorgio L. Papadopoulos.
