require(reshape2) # cat tabulate to matrix
library(misha)
library(shaman)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
source('./scripts_clean/auxFunctions.R')
options("scipen"=999)
colors <- rainbow(10)

mDBloc <- '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

obsData <- list(
        larvae_DWT     = c("hic.hic_larvae_DWT_Rep1_dm6_BeS.hic_larvae_DWT_Rep1_dm6_BeS_1",
                           "hic.hic_larvae_DWT_Rep2_dm6_BeS.hic_larvae_DWT_Rep2_dm6_BeS_1",
                           "hic.hic_larvae_DWT_Rep3_dm6_BeS.hic_larvae_DWT_Rep3_dm6_BeS_1"),
        larvae_double  = c("hic.hic_larvae_double_Rep1_dm6_BeS.hic_larvae_double_Rep1_dm6_BeS_1",
                           "hic.hic_larvae_double_Rep2_dm6_BeS.hic_larvae_double_Rep2_dm6_BeS_1",
                           "hic.hic_larvae_double_Rep3_dm6_BeS.hic_larvae_double_Rep3_dm6_BeS_1"),
        LD_gypsy2Enh   = c("hic.hic_LD_gypsy2Enh_Rep1_dm6_BeS.hic_LD_gypsy2Enh_Rep1_dm6_BeS_1",
                           "hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1",
                           "hic.hic_LD_gypsy2Enh_Rep2_dm6_BeS.hic_LD_gypsy2Enh_Rep2_dm6_BeS_1")
                       )

refSample <- "larvae_DWT"

#obsData <- list(
#	ED_PH18 = c('microC.mC_ED_PH18_Rep1_dm6_BeS.mC_ED_PH18_Rep1_dm6_BeS',
#		    'microC.mC_ED_PH18_Rep1res_dm6_BeS.mC_ED_PH18_Rep1res_dm6_BeS',
#		    'microC.mC_ED_PH18_Rep2_dm6_BeS.mC_ED_PH18_Rep2_dm6_BeS',
#		    'microC.mC_ED_PH18_Rep2res_dm6_BeS.mC_ED_PH18_Rep2res_dm6_BeS'),
#	ED_PH29 = c('microC.mC_ED_PH29_Rep1_dm6_BeS.mC_ED_PH29_Rep1_dm6_BeS',
#		    'microC.mC_ED_PH29_Rep1res_dm6_BeS.mC_ED_PH29_Rep1res_dm6_BeS',
#		    'microC.mC_ED_PH29_Rep2_dm6_BeS.mC_ED_PH29_Rep2_dm6_BeS',
#		    'microC.mC_ED_PH29_Rep2res_dm6_BeS.mC_ED_PH29_Rep2res_dm6_BeS')
#		    )

#refSample <- "ED_PH18"

samples = names(obsData)
print(paste0("Samples ",samples))

refChrom  <- "chr2L"

DomainsFile <- "./scripts_clean/list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed"
Domains     <- read.table(DomainsFile, header=F)
colnames(Domains) <- c("chrom","start","end","domain")
print(head(Domains))

RegionsFile <- "./scripts_clean/list_of_PcG_physical_domains_from_Sexton_et_al_2012_dm6.bed"
Regions     <- read.table(RegionsFile, header=F)
colnames(Regions) <- c("chrom","start","end","domain")

totalContacts <- read.table(paste0("./scripts_clean/obsContacts_per_sample_genomewide.tab"), header=F)
colnames(totalContacts) <- c("sample","count")
print(totalContacts)

if(refChrom == "chr2L")
{
    Domains     <- Domains[Domains$chrom == refChrom,]
    Regions     <- Regions[Regions$chrom == refChrom,]
}
print(head(Domains))
print(head(Regions))

outfile <- paste0("obsContacts_within_domains_in_",refChrom,"_revisions.tab")
if(!file.exists(outfile))
{
    for(sample in samples)
    {
        print(paste0("Analysing ",sample))

        result  <- data.frame()
        results <- data.frame()

        for(g in 1:nrow(Regions))
        {
            region <- Regions[g,]    

	    chrom1      <- region$chrom
    	    start1      <- region$start
    	    end1        <- region$end
	    domain1     <- region$domain
    	    print(paste(g,"of",nrow(Regions),chrom1,start1,end1,domain1,sep=" "), quote=FALSE)

	    for(d in 1:g)
	    {
	        domain <- Domains[d,]	    
	    
	        chrom2    <- domain$chrom
	        start2    <- domain$start
	        end2      <- domain$end
	        domain2   <- domain$domain
	        if(domain1 != domain2)
	        {
	            next
	        }
	        print(paste(d,"of",g,chrom2,start2,end2,domain2,sep=" "), quote=FALSE)	    

                ##### To plot scores #####
   	        interval <- gintervals.2d(chrom1,start1,end1,chrom2,start2,end2)
	        track <- obsData[[sample]]
	        print(track)
	        for(j in 1:length(track))
                {
                    vTrack <- gvtrack.create(paste0("obs",j),track[j],"area")
                }
	        data <- gextract(paste0("obs",1:length(track)),interval,iterator=interval)

	        contacts <- sum(data[,paste0("obs",1:length(track))])
	    
	        result <- data.frame(sample,chrom1,start1,end1,domain1,chrom2,start2,end2,domain2,contacts)
	        if(nrow(results) == 0)
	        {
	            results <- result
	        } else {
	            results <- rbind(results,result)
	        }
	    
	    } # Close cycle over d
        } # Close cycle over r

        colnames(results) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","contacts")
        write.table(as.matrix(results),file=outfile,sep="\t",quote=F,row.names=F,col.names=F,append=T)

    } # Close cycle over sample
}

df_obsContacts <- read.table(outfile,header=F)
colnames(df_obsContacts) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","contacts")
print(head(df_obsContacts))

result <- data.frame()
results <- data.frame()

outfile <- paste0("obsContacts_within_domains_vs_",refSample,"_in_",refChrom,"_revisions.tab")
if(!file.exists(outfile))
{
    for(domain in unique(df_obsContacts$domain1))
    {
        print(domain)

        domainContacts <- df_obsContacts[df_obsContacts$domain1 == domain,]
	print(head(domainContacts))

        chrom1     <- unique(domainContacts$chrom1)
        start1     <- unique(domainContacts$start1)
        end1       <- unique(domainContacts$end1)
        chrom2     <- unique(domainContacts$chrom2)
        start2     <- unique(domainContacts$start2)
        end2       <- unique(domainContacts$end2)    
	#print(paste0(chrom1," ",start1," ",end1," ",chrom2," ",start2," ",end2))

        refContacts <- domainContacts[domainContacts$sample == refSample,]$contacts
	refTotal <- totalContacts[totalContacts$sample == refSample,]$count	
	print(refTotal)

        for(sample in unique(domainContacts$sample))
        {   
 
            obsContacts <- domainContacts[domainContacts$sample == sample,]$contacts
	    obsTotal    <- totalContacts[totalContacts$sample == sample,]$count
	    print(obsTotal)	

	    result <- data.frame(sample,chrom1,start1,end1,domain,chrom2,start2,end2,domain,refContacts,obsContacts,log(obsContacts/refContacts)/log(2),refTotal,obsTotal,log((obsContacts/obsTotal)/(refContacts/refTotal))/log(2))

	    if(nrow(results) == 0)
            {
	        results <- result
	    } else {
	        results <- rbind(results,result)
	    }
	
        }
    }

    colnames(results) <- c("sample","chrom1","start1","end1","domain1","chrom2","start2","end2","domain2","refContacts","obsContacts",paste0("Log2FC_Contacts_",refSample),"refCoverage","obsCoverage",paste0("Log2FC_obsContactsCoverageNorm_vs_refContactsCoverageNorm_",refSample))
    write.table(as.matrix(results),file=outfile,sep="\t",quote=F,row.names=F)
}
