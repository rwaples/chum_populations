all: batch_10

batch_10: filter statistics plots  

filter: missingness hwe 

missingness: ./data/batch_10/batch_10_raw.ped ./data/batch_10/batch_10_raw.map
	plink --file ./data/batch_10/batch_10_raw --geno 0.2 --maf .1 --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_filter_loci
	plink --bfile ./data/batch_10/batch_10_filter_loci --family --geno 0.25 --mind 0.25 --maf .05 --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_filter_missingness

hwe:
	plink --bfile ./data/batch_10/batch_10_filter_missingness --write-snplist --family --hwe .0000001 midp --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_clean


statistics: basic_stats f_stats pca ibs_dist mds

basic_stats:
	plink --bfile ./data/batch_10/batch_10_clean --hardy midp --freq --family --allow-extra-chr --out ./data/batch_10/batch_10

f_stats: ./data/batch_10/batch_10_clean.bed ./data/batch_10/batch_10_clean.bim
	plink --bfile ./data/batch_10/batch_10_clean --het small-sample --ibc --fst --family --allow-extra-chr --out ./data/batch_10/batch_10

pca: 
	plink --bfile ./data/batch_10/batch_10_clean --pca header tabs var-wts --make-rel --allow-extra-chr --out ./data/batch_10/batch_10

ibs_dist: ./data/batch_10/batch_10_clean.bed ./data/batch_10/batch_10_clean.bim
	plink --bfile ./data/batch_10/batch_10_clean --pca --distance square flat-missing 1-ibs --allow-extra-chr --out ./data/batch_10/batch_10

mds: 
	Rscript ./scripts/12_mds.R ./POPINFO.txt ./data/batch_10/batch_10.mdist.id ./data/batch_10/batch_10.mdist 20 ./data/batch_10/batch_10.dist_mds


plots: plot_mds

plot_mds: 
	Rscript ./scripts/13_plot_mds.R ./data/batch_10/batch_10.dist_mds



