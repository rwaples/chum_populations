library(CCA)
library(vegan)





#rel_mat_1 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/paralogs.dom.relmatrix', sep = "", header = FALSE)
#rel_mat_2 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/on_map.codom.subsample.relmatrix', sep = "", header = FALSE)

evec_1 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/complete.codom.evec', 
                       sep = "", header = FALSE, skip =1)
evec_2 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/complete.dom.evec', 
                    sep = "", header = FALSE, skip =1)

evec_3 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/paralogs.dom.evec', 
										sep = "", header = FALSE, skip =1)

# CCA
my_cca = cc(evec_1[2:5], evec_2[2:6])


# procrustes
proc1 = procrustes( evec_2[2:7], evec_1[2:7], symmetric=TRUE)
plot(proc1, kind = 2)
plot(proc1, kind = 1, choices=c(1,2))
proc1$scale
(residuals(proc1))
proc1


proc2 = procrustes(  evec_1[2:3], evec_2[2:3], , symmetric=TRUE)
plot(proc2)
proc2$scale
proc2

proc3 = procrustes(evec_3[2:11], evec_1[2:11], symmetric=TRUE)
plot(proc3)
proc3$scale
proc3




protest(evec_2[2:6], evec_1[2:6], perm = 100000)

protest(evec_3[2:11], evec_1[2:11], perm = 100000)
					 
					 
					 
my_cca = cc(evec_1[2:6], evec_2[2:6])

my_cca$cor


plt.cc(my_cca, var.label=TRUE)

?cc
?read.table
