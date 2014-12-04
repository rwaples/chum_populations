# convert colum 6 in .ped files to 1
# convert chrom in .map file from un to 11

# run smart pca looking for informative missingness
/home/ipseg/Programs/EIGENSOFT/EIG6.0beta/bin/smartpca -p /home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.pca.missingness.par
# calculate Tracy-Widom statistics (Patterson et al. 2006)
/home/ipseg/Programs/EIGENSOFT/EIG6.0beta/bin/twstats -t /home/ipseg/Programs/EIGENSOFT/EIG6.0beta/POPGEN/twtable \
-i /home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.missingness.eval \
-o FISH_560.missingness.tw

# run smart pca looking for population structure
/home/ipseg/Programs/EIGENSOFT/EIG6.0beta/bin/smartpca -p /home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.pca.genotypes.par
# calculate Tracy-Widom statistics (Patterson et al. 2006)
/home/ipseg/Programs/EIGENSOFT/EIG6.0beta/bin/twstats -t /home/ipseg/Programs/EIGENSOFT/EIG6.0beta/POPGEN/twtable \
-i /home/ipseg/Desktop/waples/chum_populations/data/test_EIGENSTRAT/FISH_560.genotypes.eval \
-o FISH_560.genotypes.tw

