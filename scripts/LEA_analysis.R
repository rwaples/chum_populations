library("LEA")
setwd("~/Desktop/waples/chum_populations/results/batch_4/LFMM")

ped_genotypes = '/home/ipseg/Desktop/waples/chum_populations/data/batch_4/pop_genotypes/non_paralogs.ped'
lfmm_genotypes = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/non_paralogs.ped.lfmm'
geno_genotypes = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/non_paralogs.ped.geno'


# convert a .ped  file
ped2lfmm(input.file = ped_genotypes,  output.file = lfmm_genotypes)
# convert lfmm file to geno
lfmm2geno(input.file = lfmm_genotypes, geno_genotypes, force = TRUE)

# run pca
pc = pca(lfmm_genotypes, scale = TRUE, center = TRUE)
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
snmf_project = snmf(geno_genotypes, K=1:10, entropy = TRUE, repetitions = 3, CPU = 6, percentage = 0.10, project = "new")
# plot cross-entropy criterion of all runs of the project
plot(snmf_project, lwd = 5, col = "red", pch=1)


# get the cross-entropy of each run for K = 3
ce = cross.entropy(snmf_project, K = 3)
# select the run with the lowest cross-entropy
best = which.min(ce)
# We can also visualize a barplot of ancestry coeffcients as follows.
barplot(t(Q(snmf_project, K = 3, run  =best )), col = 1:8, xlab = 'Individual')

library(plyr)
pp = data.frame(Q(snmf_project, K = 3, run  =best ))
pp = arrange(pp, V4)
pp = arrange(pp, V3)
pp = arrange(pp, V2)
pp = arrange(pp, V1)
barplot(t(pp), col = 1:8)


barplot(t(order(Q(snmf_project, K = 4, run  =best ))), col = 1:4))


##LFMM##
basic_env = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/basic.env'

project = NULL
project = lfmm(lfmm_genotypes, basic_env, K = 3, repetitions = 5, CPU = 7, missing.data = TRUE, project = "new",
               iterations = 10000, burnin = 5000)
# get the zscores of each run for K = 3
zs = z.scores(project, K = 3)
# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)
# The median z-scores must be recalibrated before computing p-values. We suggest to start with computing the genomic inflation factor, λ, and then dividing the scores by λ  (Devlin and Roeder 1999)
lambda = median(zs^2)/.456
# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)

# Control of false discoveries
# To correct for multiple testing, we apply the Benjamini-Hochberg algorithm. This can be done as follows.
L = 12399
q = 0.05
w = which(sort(cp.values) < q * (1:L)/L)
candidates = order(cp.values)[w]
candidates

which.max(mlog10p.values(project, K=3, d=1, run = 2))

sort((mlog10p.values (project, K=3, d=1)))
plot(sort((mlog10p.values(project, K=3, d=1, run = 1))))

hist(cp.values)
mlog10p.values(project, K=3, d=1, run = 2)

# Put together a data frame of results
lfmm_results = data.frame('locnum' = seq(12399), 'cpvals' = cp.values, 'zscore' = zs.stouffer)
write.table(lfmm_results, file = '/home/ipseg/Desktop/waples/chum_populations/results/batch_4/LFMM/lfmm_results.txt', 
            sep ='\t', quote = FALSE, row.names = FALSE)




###
for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("expected FDR:", alpha))
  L = length(cp.values)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp.values) < alpha * (1:L) / L)
  candidates = order(cp.values)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


candidates
# best candidate is 2905 or 33643_47



###################################################
### code chunk number 7: LEA.Rnw:142-150
###################################################
# main options, K: (the number of ancestral populations), 
#        entropy: calculate the cross-entropy criterion, 
#        CPU: the number of CPUs.

# Runs with K between 1 and 10 with cross-entropy and 10 repetitions.
project = NULL
project = snmf("genotypes.geno", K=1:10, entropy = TRUE, repetitions = 10,
               project = "new")


###################################################
### code chunk number 8: LEA.Rnw:157-159
###################################################
# plot cross-entropy criterion of all runs of the project
plot(project, lwd = 5, col = "red", pch=1)


###################################################
### code chunk number 9: LEA.Rnw:165-170
###################################################
# get the cross-entropy of each run for K = 4
ce = cross.entropy(project, K = 3)

# select the run with the lowest cross-entropy
best = which.min(ce)

# get the zscores of each run for K = 3
zs = z.scores(project, K = 3)

# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)

lambda = median(zs^2)/.456


###################################################
### code chunk number 13: LEA.Rnw:229-231
###################################################
# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)


###################################################
### code chunk number 10: LEA.Rnw:189-199
###################################################
# main options, K: (the number of latent factors), 
#           CPU: the number of CPUs.

# Runs with K = 6 and 5 repetitions.
# The runs are composed of 6000 iterations including 3000 iterations
# for burnin.
# around 20 seconds per run.
project = NULL
project = lfmm("genotypes.lfmm", "gradients.env", K = 6, repetitions = 5, 
               project = "new")


###################################################
### code chunk number 11: LEA.Rnw:215-220
###################################################
# get the zscores of each run for K = 6
zs = z.scores(project, K = 6)

# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)


###################################################
### code chunk number 12: LEA.Rnw:224-225
###################################################
lambda = median(zs^2)/.456


###################################################
### code chunk number 13: LEA.Rnw:229-231
###################################################
# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)


###################################################
### code chunk number 14: LEA.Rnw:239-252
###################################################
for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("expected FDR:", alpha))
  L = length(cp.values)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp.values) < alpha * (1:L) / L)
  candidates = order(cp.values)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


###################################################
### code chunk number 15: LEA.Rnw:256-259
###################################################
# Copy of the pdf figures in the previous directory 
# for the creation of the vignette.
file.copy(list.files(".", pattern = ".pdf"), "..")

