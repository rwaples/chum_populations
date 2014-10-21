all: batch_10

batch_10: filter hwe fst ibs_matrix mds plot_mds


filter: missingness hwe 
missingness: ./data/batch_10/batch_10_raw.ped ./data/batch_10/batch_10_raw.map
	plink --file ./data/batch_10/batch_10_raw --geno 0.2 --maf .1 --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_filter_loci
	plink --bfile ./data/batch_10/batch_10_filter_loci --family --geno 0.2 --mind 0.2 --maf .05 --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_filter_missingness

hwe:
	plink --bfile ./data/batch_10/batch_10_filter_missingness --family --hwe .0000001 midp --allow-extra-chr --make-bed --out ./data/batch_10/batch_10_clean

	
fst: ./data/batch_10/batch_10_clean.bed ./data/batch_10/batch_10_clean.bim
	plink --bfile ./data/batch_10/batch_10_clean --fst --freq --family --allow-extra-chr --out ./data/batch_10/batch_10

ibs_matrix: ./data/batch_10/batch_10_clean.bed ./data/batch_10/batch_10_clean.bim
	plink --bfile ./data/batch_10/batch_10_clean --distance square flat-missing 1-ibs --allow-extra-chr --out ./data/batch_10/batch_10

mds: 
	Rscript file ./scripts/12_mds.R ./data/batch_10/batch_10.mdist
plot_mds: 
 
batch_10.pdf:

