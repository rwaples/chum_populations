import vcf 


#full_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_bt2.vcf"
#dup_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_dup_loci.vcf"

# WARNING .vcf file is 3.0 GB
vcf_file = "/media/Shared/Data/chum/populations/fb/pos_1.raw.vcf"
contigs_file = "/media/Shared/Data/chum/populations/fb/example_contigs_with_Ns"


# print basic info to screen
"vcftools --vcf {}".format(vcf_file)

# print allele frequencies to file
freq_file = "/home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/pos1.raw.freq"
"vcftools --vcf {} --freq --out {}".format(vcf_file, freq_file)

# depth per individual
ind_depth_file = "/home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/pos1.raw.ind_depth"
"vcftools --vcf {} --depth -c > {}".format(vcf_file, ind_depth_file)

# selecting variants
--mac # at least this many minor alleles called
--max-missing .5

selected_file = "/home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/pos1.selected"
"vcftools --vcf {} --max-missing 0.75 --mac 5 --minQ 200 --minDP 2 --recode --recode-INFO-all --out {}".format(vcf_file, selected_file)

# convert to BCF
"vcftools --vcf {} --contigs {} --recode-bcf --recode-INFO-all --out {}".format(vcf_file, contigs_file, vcf_file.replace(".vcf",".bcf") )




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


