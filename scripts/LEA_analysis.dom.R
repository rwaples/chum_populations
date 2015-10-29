# LFMM including duplicated loci (dominance coded)

library("LEA")
setwd("~/Desktop/waples/chum_populations/results/batch_4/LFMM")

raw_genotypes = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/complete.dom.geno'

lfmm_genotypes = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/complete.dom.lfmm'
geno_genotypes = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/complete.dom.geno'


# convert raw file to lfmm
geno2lfmm(input.file = raw_genotypes, lfmm_genotypes, force = TRUE)

# convert lfmm file to genonmf

lfmm2geno(input.file = lfmm_genotypes, geno_genotypes, force = TRUE)

# run pca
pc = pca(lfmm_genotypes, scale = FALSE, center = TRUE)
par(mfrow=c(1,1))

# Plot eigenvalues.
plot(pc, lwd=5, col="red",xlab=("PCs"),ylab="eigen")
# Plot standard deviations.
plot(pc$sdev)
# pca projection
plot(pc$projections[,1], pc$projections[,2])

# Perfom Tracy-Widom tests on all eigenvalues.
tw = tracy.widom(pc)
# display the p-values for the Tracy-Widom tests and plot the percentage of variance explained by each component
plot( -log10(tw$pvalues),tw$percentage)


# run snmf (sparse non-negative matrix factorization)
snmf_project = NULL
snmf_project = snmf(geno_genotypes, K=2:6, entropy = TRUE, repetitions = 5, CPU = 6, percentage = 0.10, project = "new", ploidy = 1)

# we may want to deal with the 

# plot cross-entropy criterion of all runs of the project
plot(snmf_project, lwd = 5, col = "red", pch=1)

# get the cross-entropy of each run for K = 3
ce = cross.entropy(snmf_project, K = 3)
# select the run with the lowest cross-entropy
best = which.min(ce)
# We can also visualize a barplot of ancestry coeffcients as follows.
barplot(t(Q(snmf_project, K = 3, run  =best )), col = 1:8, xlab = 'Individual')



##LFMM##
basic_env = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/basic.env'

lfmm_project = NULL
lfmm_project = lfmm(lfmm_genotypes, basic_env, K = 3, repetitions = 3, CPU = 7, missing.data = TRUE, project = "new", iterations = 10000, burnin = 5000)

# get the zscores of each run for K = 3
zs = z.scores(lfmm_project, K = 3)
# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)
# The median z-scores must be recalibrated before computing p-values. We suggest to start with computing the genomic inflation factor, λ, and then dividing the scores by λ  (Devlin and Roeder 1999)
lambda = median(zs^2)/.456
# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)

# Control of false discoveries
# To correct for multiple testing, we apply the Benjamini-Hochberg algorithm. This can be done as follows.
L = 35579
q = 0.05
w = which(sort(cp.values) < q * (1:L)/L)
candidates = order(cp.values)[w]
candidates

which.max(mlog10p.values(lfmm_project, K=3, d=1, run = 3))

sort((mlog10p.values(lfmm_project, K=3, d=1)))[1:100]
plot(sort((mlog10p.values(lfmm_project, K=3, d=1, run = 1))))

hist(cp.values)
mlog10p.values(lfmm_project, K=3, d=1, run = 2)

# Put together a data frame of results
lfmm_results = data.frame('locnum' = seq(35579), 'cpvals' = cp.values, 'zscore' = zs.stouffer)
write.table(lfmm_results, file = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/lfmm_results.dom.txt', sep ='\t', quote = FALSE, row.names = FALSE)
