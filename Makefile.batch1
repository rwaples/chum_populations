HWE_POP_THRESHOLD = 5
#MAF_THRESHOLD = 0.05


all: FISH_560

FISH_560: filter statistics plots  

filter: missingness hwe maf ld final

missingness:
	# remove loci not genotyped in at least 75% of indidivudals (across all pops)
	plink --file ./data/batch_1/batch_1.plink.raw --geno 0.25 --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter01
	# remove loci not genotyped in at least 75% of indidivudals (within each pop)
	plink --bfile  ./data/batch_1/batch_1.filter01 --family --geno 0.25 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter02
	# remove invididuals not genotyed at least 75% of loci	 
	plink --bfile  ./data/batch_1/batch_1.filter02 --family --mind 0.25 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter03
	# remove population 8 (CMSQUA10)
	plink --bfile  ./data/batch_1/batch_1.filter03 --family --remove-cluster-names 8 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter04
	
	#split files by family
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 1 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_01
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 2 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_02
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 3 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_03	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 4 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_04	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 5 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_05	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 6 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_06	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 7 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_07	
	#plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 8 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_08	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 9 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_09	
	plink --bfile  ./data/batch_1/batch_1.filter04 --family --keep-cluster-names 10 -allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.fam_10	

hwe:	
	# list of loci passing HWE for each pop
	plink --bfile  ./data/batch_1/batch_1.fam_01 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_01.hwe	
	plink --bfile  ./data/batch_1/batch_1.fam_02 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_02.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_03 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_03.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_04 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_04.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_05 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_05.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_06 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_06.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_07 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_07.hwe
	#plink --bfile  ./data/batch_1/batch_1.fam_08 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_08.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_09 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_09.hwe
	plink --bfile  ./data/batch_1/batch_1.fam_10 --hwe .0001 midp --hardy midp --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_10.hwe

maf:	
	# loci with MAF > .05 for each population
	plink --bfile  ./data/batch_1/batch_1.fam_01 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_01.maf	
	plink --bfile  ./data/batch_1/batch_1.fam_02 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_02.maf
	plink --bfile  ./data/batch_1/batch_1.fam_03 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_03.maf
	plink --bfile  ./data/batch_1/batch_1.fam_04 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_04.maf
	plink --bfile  ./data/batch_1/batch_1.fam_05 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_05.maf
	plink --bfile  ./data/batch_1/batch_1.fam_06 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_06.maf
	plink --bfile  ./data/batch_1/batch_1.fam_07 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_07.maf
	#plink --bfile  ./data/batch_1/batch_1.fam_08 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_08.maf
	plink --bfile  ./data/batch_1/batch_1.fam_09 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_09.maf
	plink --bfile  ./data/batch_1/batch_1.fam_10 --maf .05 --allow-extra-chr --allow-no-sex --write-snplist --out ./data/batch_1/batch_1.fam_10.maf
	
	# list of loci passing HWE in at least $(HWE_POP_THRESHOLD) populations
	python ./scripts/11a_plink_HWE.py $(HWE_POP_THRESHOLD) ./data/batch_1/batch_1.fam_01.hwe.snplist ./data/batch_1/batch_1.fam_02.hwe.snplist ./data/batch_1/batch_1.fam_03.hwe.snplist ./data/batch_1/batch_1.fam_04.hwe.snplist ./data/batch_1/batch_1.fam_05.hwe.snplist ./data/batch_1/batch_1.fam_06.hwe.snplist ./data/batch_1/batch_1.fam_07.hwe.snplist ./data/batch_1/batch_1.fam_09.hwe.snplist ./data/batch_1/batch_1.fam_10.hwe.snplist
	# list of loci passing MAF threshold in at least one population	
	python ./scripts/11b_plink_MAF.py ./data/batch_1/batch_1.fam_01.maf.snplist ./data/batch_1/batch_1.fam_02.maf.snplist ./data/batch_1/batch_1.fam_03.maf.snplist ./data/batch_1/batch_1.fam_04.maf.snplist ./data/batch_1/batch_1.fam_05.maf.snplist ./data/batch_1/batch_1.fam_06.maf.snplist ./data/batch_1/batch_1.fam_07.maf.snplist ./data/batch_1/batch_1.fam_09.maf.snplist ./data/batch_1/batch_1.fam_10.maf.snplist

	# remove loci failing the HWE or MAF filters
	plink --bfile  ./data/batch_1/batch_1.filter04 --extract ./data/FISH_560/passing_HWE.snps --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter05
	plink --bfile  ./data/batch_1/batch_1.filter05 --extract ./data/FISH_560/passing_MAF.snps --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter06
	
ld:
	# LD pruning (LD stats don't work per family)
	#plink --bfile ./data/batch_1/batch_1.filter06 --family --ld-window 5 --r2 -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.filter06
	plink --bfile ./data/batch_1/batch_1.filter06 --family --indep-pairwise 3 1 .2 -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.filter06.pruned
	plink --bfile  ./data/batch_1/batch_1.filter06 --extract ./data/batch_1/batch_1.filter06.pruned.prune.in --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter07

final:
	# try a pruning step
	plink --bfile  ./data/batch_1/batch_1.filter07 --bp-space 70 --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.filter08	
	plink --bfile  ./data/batch_1/batch_1.filter08 --allow-extra-chr --allow-no-sex --make-bed --out ./data/batch_1/batch_1.final
	plink --bfile  ./data/batch_1/batch_1.final --allow-extra-chr --allow-no-sex --recode 12 A --out ./data/batch_1/batch_1.final


statistics: basic_stats f_stats count_retained pca ibs_dist mds

basic_stats:
	plink --bfile ./data/batch_1/batch_1.final --family --freq --missing -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.final

f_stats:
	plink --bfile ./data/batch_1/batch_1.final --family --het small-sample --ibc --fst -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.final
	plink --bfile ./data/batch_1/batch_1.final --family --fst -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.final

# multivariate stats
count_retained:
	wc -l < ./data/batch_1/batch_1.final.bim > ./data/batch_1/batch_1.num_loci_retained
	wc -l < ./data/batch_1/batch_1.final.fam > ./data/batch_1/batch_1.num_inds_retained

LOCI_RETAINED := $(shell cat ./data/batch_1/batch_1.num_loci_retained)
INDS_RETAINED := $(shell cat ./data/batch_1/batch_1.num_inds_retained)

pca:
	plink --bfile ./data/batch_1/batch_1.final --pca $(INDS_RETAINED) header tabs var-wts --make-rel -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.final

ibs_dist: 
	plink --bfile ./data/batch_1/batch_1.final --distance square flat-missing 1-ibs -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.final

mds:
	# two different ways of running an MDS
	Rscript ./scripts/12_mds.R ./data/POPINFO.txt ./data/batch_1/batch_1.final.mdist.id ./data/batch_1/batch_1.final.mdist 20 ./data/batch_1/batch_1.dist_mds
	#plink --bfile ./data/batch_1/batch_1.filter06 --cluster --mds-plot 3  eigvals -allow-extra-chr --allow-no-sex --out ./data/batch_1/batch_1.filter06
	
plots: plot_mds

plot_mds: 
	Rscript ./scripts/13_plot_mds.R ./data/batch_1/batch_1.dist_mds



