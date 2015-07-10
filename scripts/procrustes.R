library(CCA)
library(vegan)


rel_mat_1 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/paralogs.dom.relmatrix', 
    sep = "", header = FALSE)
rel_mat_2 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/on_map.codom.subsample.relmatrix', sep = "", header = FALSE)

evec_1 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/on_map.codom.subsample.evec', 
                       sep = "", header = FALSE, skip =1)
evec_2 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/on_map.codom.evec', 
                    sep = "", header = FALSE, skip =1)

# CCA
my_cca = cc(evec_1[2:5], evec_2[2:6])


# procrustes
proc1 = procrustes( evec_2[2:5], evec_1[2:5], symmetric=TRUE)
plot(proc1, kind = 1, choices=c(1,2))
proc1$scale
hist(residuals(proc1))
proc1


proc2 = procrustes(  evec_1[2:3], evec_2[2:3], , symmetric=TRUE)
plot(proc2)
proc2$scale
proc2

protest(evec_2[2:6], evec_1[2:6], perm = 100000)


my_cca = cc(evec_1[2:6], evec_2[2:6])

my_cca$cor


plt.cc(my_cca, var.label=TRUE)

?cc
?read.table
