import vcf 


full_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_bt2.vcf"
dup_vcf_file = "/home/ipseg/Desktop/waples/chum_populations/batch_02_dup_loci.vcf"

# list loads the whole thing in memory, not always necessary
dup_vcf = list(vcf.Reader(open(dup_vcf_file)))

# histogram of number of call genotypes
hist([rec.num_called for rec in dup_vcf])


# histogram of depth (DP) per locus
hist([rec.INFO['DP'] for rec in dup_vcf], bins = 100)

# dont do this, as ABP needs to be flattened first
hist([rec.INFO['ABP'] for rec in dup_vcf], bins = 100)

dup_vcf[0].INFO['ABP']