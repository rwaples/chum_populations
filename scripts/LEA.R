# get bioconductor

source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
browseVignettes("LEA")

# take from :  http://www.bioconductor.org/packages/release/bioc/vignettes/LEA/inst/doc/LEA.R

### R code from vignette source 'LEA.Rnw'

###################################################
### code chunk number 1: LEA.Rnw:48-53
###################################################
# creation of a directory for the LEA analyses


# set this new directory as the working directory
setwd("~/Desktop/waples/chum_populations/results/LEA_analyses")


###################################################
### code chunk number 2: LEA.Rnw:77-88
###################################################
library(LEA)
# Creation of the genotypic file "genotypes.lfmm"
# with 400 SNPs for 50 individuals.
data("tutorial")
# in the lfmm format
write.lfmm(tutorial.R, "genotypes.lfmm")
# in the geno format
write.geno(tutorial.R, "genotypes.geno")
# creation of the environment file, gradient.env.
# It contains 1 environmental variable for 50 individuals.
write.env(tutorial.C, "gradients.env")


###################################################
### code chunk number 3: LEA.Rnw:99-110
###################################################
# run of pca
# Available options, K (the number of PCs calculated), 
#                    center and scale. 
# Creation of   genotypes.pcaProject - the pcaProject object.
#               a directory genotypes.pca containing:
# Create files: genotypes.eigenvalues - eigenvalues,        
#               genotypes.eigenvectors - eigenvectors,
#               genotypes.sdev - standard deviations,
#               genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("genotypes.lfmm", scale = TRUE)


###################################################
### code chunk number 4: LEA.Rnw:115-118
###################################################
# Perfom Tracy-Widom tests on all eigenvalues.
# create file: tuto.tracyWidom - tracy-widom test information.  
tw = tracy.widom(pc)


###################################################
### code chunk number 5: LEA.Rnw:121-123
###################################################
# display the p-values for the Tracy-Widom tests. 
tw$pvalues[1:10]


###################################################
### code chunk number 6: LEA.Rnw:128-130
###################################################
# plot the percentage of variance explained by each component
plot(tw$percentage)


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
ce = cross.entropy(project, K = 4)

# select the run with the lowest cross-entropy
best = which.min(ce)


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

