# TODO

# use the same set of indidivudals for each analysis


library(nFactors)
library(stringr)
# need to remove A1\t from the eigenvec.var file
locus_weights <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/FISH_560/FISH_560.final.eigenvec.var", header = TRUE)

eigenvec <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/FISH_560/FISH_560.final.eigenvec", sep = "\t", header = TRUE)
eigenval <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/FISH_560/FISH_560.final.eigenval", sep = "\t", header = FALSE)

smartPCA_missingness_eval <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.missingness.eval", sep = "\t", header = FALSE)
smartPCA_missingness_evec <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.missingness.evec", header = FALSE, skip =  1, row.names= 1)

smartPCA_genotypes_eval <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.genotypes.eval", sep = "\t", header = FALSE)
smartPCA_genotypes_evec <- read.table("/home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.genotypes.evec", header = FALSE, skip =  1, row.names= 1)

names(smartPCA_genotypes_evec) <- paste('dim', seq(1, length(smartPCA_genotypes_evec)), sep = "")
names(smartPCA_missingness_evec) <-  paste('dim', seq(1, length(smartPCA_missingness_evec)), sep = "")
missingness_IDs = data.frame(str_split_fixed(row.names(smartPCA_missingness_evec), ":", 2))
genotype_IDs = data.frame(str_split_fixed(row.names(smartPCA_genotypes_evec), ":", 2))
names(missingness_IDs) <- c('pop', 'ind')
names(genotype_IDs) <- c('pop', 'ind')
smartPCA_missingness <- cbind(missingness_IDs, smartPCA_missingness_evec)
smartPCA_genotypes <- cbind(genotype_IDs, smartPCA_genotypes_evec)

smartPCA_genotypes_eval$pc_axis <- seq(nrow(smartPCA_genotypes_eval))

# simulate a parallele analysis, using: nFactors.parallel()
PA_sim <- parallel(subject = nrow(eigenvec), var = length(eigenvec)-2, rep = 1000, quantile = .01)
PA_sim <- parallel(subject = nrow(eigenvec), var = length(smartPCA_evec)-1, rep = 1000, quantile = .01)
#PA_sim$eigen$mevpea = mean of eigenvalues distribution
#PA_sim$eigen$sevpea = sd of eigenvalues distribution
#PA_sim$eigen$qevpea = quantile of eigenvalues distribution
#PA_sim$eigen$sqevpea = std err of the quantile of eigenvalues distribution

plotParallel(PA_sim, x= eigenval$V1)
plotParallel(PA_sim, x= smartPCA_genotypes_eval$V1)
plotuScree(x= eigenval$V1)

library(ggplot2)
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}

ggplot(data = smartPCA_missingness, aes(x = dim1, y = dim2, color = pop)) + geom_point(size = 6, alpha = .8) + coord_equal() +theme_bw() + 
  xlab(paste("dim1: ", smartPCA_missingness_eval[,1][1], "%") )  + ylab(paste("dim2: ", smartPCA_missingness_eval[,1][2], "%") ) 
ggplot(data = smartPCA_missingness, aes(x = dim1, y = dim3, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()

ggplot(data = smartPCA_genotypes, aes(x = dim1, y = dim2, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()  + 
  xlab(paste("dim1: ", smartPCA_genotypes_eval[,1][1], "%") )  + ylab(paste("dim2: ", smartPCA_genotypes_eval[,1][2], "%") ) 
ggplot(data = smartPCA_genotypes, aes(x = dim3, y = dim4, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()
ggplot(data = smartPCA_genotypes, aes(x = dim5, y = dim6, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()
ggplot(data = smartPCA_genotypes, aes(x = dim7, y = dim8, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()
ggplot(data = smartPCA_genotypes, aes(x = dim9, y = dim10, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()
ggplot(data = smartPCA_genotypes, aes(x = dim11, y = dim12, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()
ggplot(data = smartPCA_genotypes, aes(x = dim2, y = dim3, color = pop)) + geom_point(size = 6, alpha = .8) + theme_bw()

#scree plot
ggplot(data = smartPCA_genotypes_eval, aes(y=V1, x=pc_axis )) + geom_point(shape = 16, size = 4, alpha = .8) + geom_path(size = 1, alpha = .2) + theme_bw()

cor(smartPCA_missingness$dim4[1:100], smartPCA_genotypes$dim4)
for (xx in seq(3, 100)) {
  print(cor(smartPCA_missingness[xx][1:158,], smartPCA_genotypes[xx]))
}


names(smartPCA_genotypes_evec)
plot(smartPCA_evec[,1], smartPCA_evec[,2])
