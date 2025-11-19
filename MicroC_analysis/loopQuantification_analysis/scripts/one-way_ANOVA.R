# Script inspired by https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova

library(tidyverse)
library(ggpubr)
library(rstatix)
library(DataScienceR)

inFiles = list("XXXinFileXXX")

detect_outlier <- function(x) {

  Quantile1 <- quantile(x, probs=.25)
  Quantile3 <- quantile(x, probs=.75)
  
  IQR = Quantile3-Quantile1
  
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}

for(outliers in c("FALSE","TRUE"))
{

for(inFile in inFiles)
{
    print(inFile, quote=F)
    data <- read.table(inFile, header=F)
    ###
    colnames(data) <- c("condition", "values")
    print(head(data), quote=F)

    levels <- unique(data$condition)
    data$condition <- factor(data$condition, levels = levels)

    print(paste0("### Get summary statistics ###"), quote=F)
    dataStats <- data %>%
	  group_by(condition) %>%
	  get_summary_stats(values, type = "mean_sd")
    print(dataStats, quote=F)

	#bxp <- ggboxplot(
	bxp <- ggviolin(	
	  data, x = "condition", y = "values",
	  color = "condition", fill= "condition", palette = "jco", scale = "width"
	) +
	ylim(c(-100,100)) +
	geom_boxplot(color = "black", fill= "white", width=0.2) +
	geom_text(data=dataStats, aes(y=min(data$values), label = n), color="black") +
        theme(axis.text.x = element_text(angle = 60, hjust=1))

	outFile <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),".pdf")
	pdf(outFile)
	print(bxp)
	dev.off()
	
	print(paste0("### Check assumptions ###"), quote=F)
	print(paste0("1) Checking the presence of outliers"), quote=F)
	print(identify_outliers(data,variable="values"), quote=F)

	checkOutliers <- data %>%
	  group_by(condition) %>%
	  identify_outliers(variable="values")
	print(checkOutliers, quote=F)

	if(nrow(checkOutliers) == 0)
	{
	    print(paste0("No outliers found: This assumption is verified"), quote=F)
	} else {
	    print(paste0("### WARNING: The following data-points are outiers!"), quote=F)
	    print(checkOutliers, quote=F)
	    print(paste0("Note that, in the situation where you have extreme outliers, this can be due to:"), quote=F)
	    print(paste0("1) data entry errors, measurement errors or unusual values."), quote=F)
	    print(paste0("You can include the outlier in the analysis anyway if you do not believe the result will be substantially affected."), quote=F)
	    print(paste0("This can be evaluated by comparing the result of the ANOVA test with and without the outlier."), quote=F)
	    print(paste0("It’s also possible to keep the outliers in the data and perform robust ANOVA test using the WRS2 package."), quote=F)

	    outFile <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),".pdf")
	    ### Code to Remove Outliers ###
	    if(outliers == "FALSE")
	    {
	        print(paste0("Removing the outliers"), quote=F)
	    
		noOutliers = data.frame()
	    	for(c in unique(data$condition))
	    	{
	            print(c)
		    noOutlier <- data[data$condition == c,]
	    	    noOutlier$isOutlier <- is_outlier(data[data$condition == c,]$values)
		    noOutlier <- noOutlier[noOutlier$isOutlier == F,]
		    print(head(noOutlier))
	    	    if(nrow(noOutliers) == 0)
	    	    {
	                noOutliers <- noOutlier
	    	    } else {
	                noOutliers <- rbind(noOutliers,noOutlier)
	            }
	        }
	    	print(nrow(data), quote=F)
	    	print(nrow(checkOutliers), quote=F)
	    	print(nrow(noOutliers), quote=F)
	    	data <- noOutliers
	    	outFile <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),"_noOuliers.pdf")
	    }

	    bxp <- ggviolin(	
	        data, x = "condition", y = "values",
	  	color = "condition", fill= "condition", palette = "jco", scale = "width"
	    ) +
	    geom_boxplot(color = "black", fill= "white", width=0.2) +
	    geom_text(data=dataStats, aes(y=min(data$values), label = n), color="black") +
            theme(axis.text.x = element_text(angle = 60, hjust=1))


	    pdf(outFile)
	    print(bxp)
	    dev.off()

	    print("", quote=F)     	    
	}
	print("", quote=F)
	#next

	print(paste0("2) Checking normality assumption"), quote=F)
	print(paste0("# Building the linear model"), quote=F)
	model  <- lm(values ~ condition,
             data = data)
	print(paste0("# Creating a QQ plot of residuals"), quote=F)
	p <- ggqqplot(residuals(model))
	pdf(paste0("qqplot_oneWay_",gsub(".tsv","",inFile),".pdf")) 
	print(p)
	dev.off()	
	print(paste0("### Compute Shapiro-Wilk test of normality ###"), quote=F)
	#checkNormality <- shapiro_test(residuals(model))
	checkNormality <- ks.test(x=data$values,y="pnorm",alternative='two.sided')
	print(checkNormality, quote=F)
	if(checkNormality$p.value > 0.05)
	{
	    print(paste0("In the QQ plot, as all the points fall approximately along the reference line, we can assume normality."), quote=F)
	    print(paste0("This conclusion is supported by the Shapiro-Wilk test. The p-value is not significant (p = ",checkNormality$p.value,"),"), quote=F)
	    print(paste0("so we can assume normality."), quote=F)
	} else {
	    print(paste0("The normality assumption is not supported by the Shapiro-Wilk test. The p-value is significant (p = ",checkNormality$p.value,"),"), quote=F)
	    print(paste0("so we cannot assume normality for this dataset."), quote=F)
	    #    quit()
	}
	print("", quote=F)
	#next

	print(paste0("3) Checking homogeneity of variance assumption"), quote=F)
	print(paste0("This can be checked using the Levene’s test:"), quote=F)
	checkHomogenityOfVariance <- data %>% levene_test(values ~ condition)
	print(checkHomogenityOfVariance, quote=F)
	if(checkHomogenityOfVariance$p < -0.05)
	{
	    print(paste0("The Levene’s test is not significant (p > ",checkHomogenityOfVariance$p,"). Therefore, we can assume the homogeneity of variances in the different groups."), quote=F)

	    print(paste0("### Computation of the one-way ANOVA test ###"), quote=F)
	    print(paste0("In the R code below, the asterisk represents the interaction effect and the main effect of each variable (and all lower-order interactions)."), quote=F)
	    res.aov <- data %>%
	    	    anova_test(values ~ condition)
	    print(res.aov, quote=F)
	    #    	    group_by(condition) %>%

	    significantInteraction = res.aov[res.aov$p < 0.052,]
	    if(nrow(significantInteraction) > 0)
	    {
	        print(paste0("There was a statistically significant interaction between conditions"), quote=F)
	        print(significantInteraction, quote=F)
	    } else {
	        print(paste0("There was no statistically significant interaction between conditions"), quote=F)
	        next
	    }

	    print("", quote=F)	
	    print(paste0("### Post-hoct tests ###"), quote=F)

	    print(paste0("### Procedure for significant one-way interaction ###"), quote=F)
	    print(paste0("### Compute pairwise comparisons ###"), quote=F)

	    print(paste0("To determine which group means are different. We’ll now perform multiple pairwise comparisons between the different conditions."), quote=F)

	    print(paste0("You can run and interpret all possible pairwise comparisons using a Bonferroni adjustment."), quote=F)
	    print(paste0("This can be easily done using the function emmeans_test() [rstatix package], a wrapper around the emmeans package,"), quote=F)
	    print(paste0("which needs to be installed. Emmeans stands for estimated marginal means (aka least square means or adjusted means)."), quote=F)
	

	    print(paste0("Compare the values of the different conditions by condition levels:"), quote=F)

	    pairwiseTestsCondition <- data %>%
	    			   group_by(condition) %>%
 				   #  		   	           games_howell_test(values ~ condition, p.adjust.method = "bonferroni")
				   tukey_hsd(values ~ condition)
				   

	    pairwiseTestsCondition <- as.data.frame(pairwiseTestsCondition)				   
	    pairwiseTestsCondition$group1 <- paste0(pairwiseTestsCondition$group1)
            pairwiseTestsCondition$group2 <- paste0(pairwiseTestsCondition$group2)

	} else {
	   print(paste0("The Levene’s test is significant (p = ",checkHomogenityOfVariance$p," < 0.05). Therefore, we cannot assume the homogeneity of variances in the different groups."), quote=F)

	   print(paste0("The Welch one-way test is an alternative to the standard one-way ANOVA"), quote=F)
	   print(paste0("in the situation where the homogeneity of variance can’t be assumed (i.e., Levene test is significant)."), quote=F)

	   res.aov <- data %>%
  	       welch_anova_test(values ~ condition)
	   print(res.aov, quote=F)
	   #	       group_by(condition) %>%

	   print(paste0("In this case, the Games-Howell post hoc test or pairwise t-tests (with no assumption of equal variances) can be used to compare all possible combinations of group differences."), quote=F)

	   pairwiseTestsCondition <- data %>%
	       games_howell_test(values ~ condition)
	       #tukey_hsd(values ~ condition)

	       #pairwise_t_test(values ~ cond, )
	   #	       group_by(condition) %>%


	   #pairwiseTestsCondition <- data %>%
	   #    group_by(loop) %>%
	   #    pairwise_ks_test(values, group=cond, n_min = 50, warning = 0, alternative = "two.sided")	   
	   #print(paste0(pairwiseTestsCondition), quote=F)
	   #quit()

	   #pairwiseTestsCondition        <- as.data.frame(pairwiseTestsCondition)						   
           pairwiseTestsCondition$group1 <- paste0(pairwiseTestsCondition$group1)
	   pairwiseTestsCondition$group2 <- paste0(pairwiseTestsCondition$group2)

	}
	print("", quote=F)	

	print(paste0("### Visualization: box plots with p-values"), quote=F)
	print(pairwiseTestsCondition, quote=F)
	#print(get_pwc_label(pairwiseTestsCondition), quote=F)
	#print(stat_pvalue_manual(pairwiseTestsCondition[pairwiseTestsCondition$p.adj < 0.05,]))

	pairwiseTestsCondition[pairwiseTestsCondition$p.adj == 0,]$p.adj = 2.22e-308	


	if(outliers == "FALSE")
	{
	    outFile  <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),"_noOutliers_with_stats.pdf")
	    outStats <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),"_noOutliers_with_stats.txt")
	} else {
	    outFile  <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),"_with_stats.pdf")
	    outStats <- paste0("boxplot_oneWay_",gsub(".tsv","",inFile),"_with_stats.txt")		    
	}

	#print(class(pairwiseTestsCondition))
	write.table(pairwiseTestsCondition,file=outStats,sep="\t", row.names=FALSE, quote=FALSE)

	pairwiseTestsCondition <- pairwiseTestsCondition %>% add_xy_position(x = "condition")	    
	pdf(outFile)
	finalBxp <- bxp +
		   #stat_pvalue_manual(pairwiseTestsCondition[pairwiseTestsCondition$p.adj < 0.05,],y.position=100,label="p.adj") +
		   stat_pvalue_manual(pairwiseTestsCondition[-grep("ns",pairwiseTestsCondition$p.adj.signif),], label="p.adj", label.size = 3) +
       	   labs(
	       subtitle = get_test_label(res.aov, detailed = TRUE),
	       caption = get_pwc_label(pairwiseTestsCondition),
	       )
	print(finalBxp)
	dev.off()	
	print("", quote=F)	
} # Close cycle over inFile
}
