import vcf 


#full_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_bt2.vcf"
#dup_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_dup_loci.vcf"
vcf_file = "/home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.raw.vcf"

# list loads the whole thing in memory, not always necessary


my_vcf = list(vcf.Reader(open(vcf_file)))
hist([rec.INFO['AC'] for rec in my_vcf])


# histogram of number of call genotypes
hist([rec.num_called for rec in my_vcf])


# histogram of depth (DP) per locus
hist([rec.INFO['DP'] for rec in dup_vcf], bins = 100)

# dont do this, as ABP needs to be flattened first
hist([rec.INFO['ABP'] for rec in dup_vcf], bins = 100)

dup_vcf[0].INFO['ABP']

# using vcftools
vcftools --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.raw.vcf \
--out /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all \
--max-missing 0.75 \
--minQ 100 \
--maf 0.05 \
--recode

vcftools --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.recode.vcf \
--out /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all \
--freq

vcftools --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.recode.vcf \
--out /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all \
--plink-tped


