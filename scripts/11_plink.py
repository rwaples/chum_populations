PLINK_EXE = "Z:/Waples/FISH_560/programs/plink/plink.exe"

Z:/Waples/FISH_560/programs/plink/plink.exe --file "Z:/Waples/FISH_560/proj_data/plink/batch_10.plink" --geno 0.2 --maf .1 --allow-extra-chr --make-bed --out batch_10_filter_loci

Z:/Waples/FISH_560/programs/plink/plink.exe --bfile "Z:/Waples/FISH_560/proj_data/plink/batch_10_filter_loci" --geno 0.2 --mind 0.2 --maf .1 --distance square flat-missing 1-ibs --allow-extra-chr --make-bed --out batch_10_filter_inds 

Z:/Waples/FISH_560/programs/plink/plink.exe --bfile "Z:/Waples/FISH_560/proj_data/plink/batch_10_filter_inds" --fst --freq --family --out batch_10_filtered --allow-extra-chr