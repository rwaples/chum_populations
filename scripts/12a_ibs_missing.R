#missingness dist matrices


raw_data <- read.table("C:/Users/Kele/Desktop/chum_populations/data/FISH_560/FISH_560.final.raw", header = TRUE)

geno_data <- data.matrix(raw_data[,7:length(raw_data)])
missing_data = geno_data
missing_data[is.na(geno_data)] <- 1
missing_data[!is.na(geno_data)] <- 0

nspns <- ncol(raw_data)- 6
ninds <- nrow(raw_data)

genotype_dist <- matrix(nrow = ninds, ncol = ninds)
missingness_dist <- matrix(nrow = ninds, ncol = ninds)

for (indx_a in seq(ninds)){
  for (indx_b in seq(ninds)){
    bb <- abs(as.numeric(geno_data[indx_a,]) - as.numeric(geno_data[indx_b,]))
    genotype_dist[indx_a, indx_b] <- mean(bb, na.rm = TRUE)/2
  }
}

for (indx_a in seq(ninds)){
  for (indx_b in seq(ninds)){
    bb <- abs(as.numeric(missing_data[indx_a,]) - as.numeric(missing_data[indx_b,]))
    missingness_dist[indx_a, indx_b] <- mean(bb, na.rm = TRUE)/2
  }
}

library(vegan)
mant <- mantel(genotype_dist, missingness_dist, method = "pearson", permutations = 10000)
#proc <- protest(genotype_dist, missingness_dist)
#plot(proc, kind = 1)

# CADM
library(ape)

genotype_dist[genotype_dist == 0] <- NA
missingness_dist[missingness_dist == 0] <- NA

mean(genotype_dist, na.rm = TRUE)
min(genotype_dist, na.rm = TRUE)
max(genotype_dist, na.rm = TRUE)  

mean(missingness_dist, na.rm = TRUE)
min(missingness_dist, na.rm = TRUE)
max(missingness_dist, na.rm = TRUE)