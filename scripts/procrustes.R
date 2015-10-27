library(CCA)
library(vegan)





#rel_mat_1 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/paralogs.dom.relmatrix', sep = "", header = FALSE)
#rel_mat_2 = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/on_map.codom.subsample.relmatrix', sep = "", header = FALSE)

evec_a = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/paralogs.dom.evec',                        sep = "", header = FALSE, skip =1)

evec_b = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/complete.dom.evec', 
                       sep = "", header = FALSE, skip =1)
evec_c = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/EIGENSOFT/complete.dom.subsample.evec', 
                    sep = "", header = FALSE, skip =1)


# CCA
my_cca = cc(evec_1[2:5], evec_2[2:6])


# procrustes
proc_ab = procrustes( evec_a[2:7], evec_b[2:7], symmetric=TRUE, )
summary(proc_ab)
plot(proc_ab, kind = 2)
plot(proc_ab, kind = 1, choices=c(1,2))

prot_ab = protest(evec_a[2:7], evec_b[2:7], permutations = how(nperm = 10000))
prot_ac = protest(evec_a[2:7], evec_c[2:7], permutations = how(nperm = 10000))
prot_bc = protest(evec_b[2:7], evec_c[2:7], permutations = how(nperm = 10000))

prot_ab
prot_ac
prot_bc


summary(prot_bc)
protest(evec_a[2:7], evec_b[2:7], permutations = how(nperm = 10000))
protest(evec_a[2:7], evec_b[2:7], permutations = how(nperm = 10000))

#protest(proc_ab[2:7], proc_ab[2:7], perm = 10000, scores = "sites")

proc_ab$scale
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
