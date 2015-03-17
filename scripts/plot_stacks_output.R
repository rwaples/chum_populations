# plot catalog
library(ggplot2)
library(plyr)


tags_file <- '/home/ipseg/Desktop/waples/chum_populations/data/batch_1/catalog/batch_1.catalog.tags.tsv'
snps_file <- '/home/ipseg/Desktop/waples/chum_populations/data/batch_1/catalog/batch_1.catalog.snps.tsv'
alleles_file <- '/home/ipseg/Desktop/waples/chum_populations/data/batch_1/catalog/batch_1.catalog.alleles.tsv'

tags <- read.table(tags_file, sep ="\t", header = FALSE)
snps <- read.table(snps_file, sep ="\t", header = FALSE)
alleles <- read.table(alleles_file, sep ="\t", header = FALSE)

names(tags) <- c("sql_ID", "batch", "locus", "ref", "pos", "strand", "seq_type", "unused", "merged", "seq", "deleveraged", "blacklisted", "lumberjack", "loglike" )
names(snps) <-  c("sql_ID", "batch", "locus", 'pos', "type", "likelihood_ratio", "nuc_1", "nuc_2", "nuc_3", "nuc_4")
names(alleles) <- c("sql_ID", "batch", "locus", "haplotype", "percent", "count")

plot_snps_per_tag <- function(alleles.df){

}
  
# alleles per variable locus
alleles.count = ddply(alleles, 'locus', summarize, allele_count = length(locus))
ggplot(data=alleles.count[alleles.count$count < 11,], aes(x=allele_count))   + geom_histogram(binwidth=1) +  xlim(c(2,10)) + scale_x_discrete(breaks = seq(2,10) )

# snps per locus
snps.count = ddply(snps, 'locus', summarize, snp_count = length(locus))
ggplot(data=snps.count[snps.count$count < 11,], aes(x=snp_count))   + geom_histogram(binwidth=1) +  xlim(c(1,10)) + scale_x_discrete(breaks = seq(1,10))

# 
tags.count = ddply(tags, 'ref', summarize, tag_count = length(ref))
ggplot(data=tags.count[tags.count$tag_count < 11,], aes(x=tag_count)) + geom_histogram( binwidth=1) +  xlim(c(1,10)) + scale_x_discrete(breaks = seq(1,10))


plot_snps_per_tag(alleles)






sumstats_summary_allsites <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/populations_output/batch_2.sumstats_summary.allsites.tsv", sep = '\t', header = TRUE)
names(sumstats_summary_allsites)

sumstats_summary_variants <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/populations_output/batch_2.sumstats_summary.variants.tsv",sep = '\t',header = TRUE)
names(sumstats_summary_variants)

ggplot(data = sumstats_summary_variants, aes(x = Pop.ID)) + geom_bar(stat="identity", aes(y=Pi))
  geom_bar(stat="identity", aes(y=Obs.Hom/Exp.Hom))

ggplot(data = sumstats_summary_allsites, aes(x = Pop.ID)) + geom_bar(stat="identity", aes(y=Pi))


sumstats <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/populations_output/batch_2.sumstats.tsv",sep = '\t', header = TRUE)
names(sumstats)
ggplot(data=sumstats, aes(x=Pop.ID)) + geom_violin(adjust = 1, alpha = .5, aes(y=Pi, fill = Pop.ID))  

ggplot(data=sumstats, aes(x=Locus.ID)) + geom_jitter(alpha = .5, aes(y=P, group = Pop.ID, color = Pop.ID))  

ibc <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.final.ibc", header = TRUE)
ggplot(data=ibc, aes(x=factor(FID))) + geom_point(alpha = .5, aes(y=Fhat1, group = FID))

ggplot(data=ibc, aes(x=NOMISS)) + geom_point(alpha = .5, aes(y=Fhat2, color = factor(FID)), size = 4)  
ggplot(data=ibc, aes(x=Fhat2)) + geom_point(alpha = .5, aes(y=Fhat3, color = factor(FID)), size = 4)  

frq.strat <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.final.frq.strat",header = TRUE)
ggplot(data=frq.strat, aes(x=factor(CLST))) + geom_violin(adjust = 1, alpha = .5, aes(y=MAF, fill = factor(CLST))) 
ggplot(data=frq.strat, aes(x=MAF)) + geom_histogram(binwidth = .03) 
ggplot(data=frq.strat, 

ggplot(data=frq.strat, aes(x=MAF, y = MAF)) + geom_point(alpha = .5, aes(group = factor(CLST))) + facet_grid(factor(CLST)~factor(CLST))

library(reshape2)
library(GGally)
MAFcast <- dcast(frq.strat, SNP~CLST, value.var= "MAF")

names(MAFcast) <- c('one', 'two', 'three','four','five','six','seven','eight','nine','ten')
ggpairs(MAFcast, 2:9, lower = list(continuous = "points"), alpha =.2)






ggplot(data=MAFcast, aes(x=1, y = 2)) + geom_point(alpha = .5, aes(group = factor(CLST))) + facet_grid(factor(CLST)~factor(CLST))

plotmatrix(MAFcast[2:9])

ggplot(data=MAFcast)


pairs(MAFcast)
ddply(frq.strat, c(SNP,CLST), mutate,, 


ggplot(data=frq.strat, aes(x=factor(CLST))) + geom_line(aes(y=MAF, group = SNP), alpha = .1) 

hwe <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.filter09.hwe",header = TRUE)
ggplot(data=hwe) + geom_point(aes(y = sort(P), x = sort(runif(length(P), max = 1)))) + geom_abline(intercept = 0)
ggplot(data=hwe, aes(x = P)) + geom_histogram(binwidth=.01) 
ggplot(data=hwe, aes(sample = P)) + geom_point(stat='qq')





punif(.2)



source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(SimRAD)

simseq <- sim.DNAseq(size=300000000, GCfreq=0.51)
#SbfI
cs_5p1 <- "CCTGCA"
cs_3p1 <- "GG"
simseq.dig <- insilico.digest(simseq, cs_5p1, cs_3p1, verbose=TRUE)

?sim.DNAseq
?insilico.digest

