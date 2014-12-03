import glob
import os
import collections
import sys

snplist_files = sys.argv[2:]
HWE_POP_THRESHOLD = int(sys.argv[1])

# build a list of all the snps from each family.  
# snps will be present x times if they appear in x families
concatenated_snps = list()
for snplist_file in snplist_files:
    with open(snplist_file, 'r') as INFILE:
        snps = [xx.strip() for xx in INFILE.readlines()]
        concatenated_snps += snps

snp_counts = collections.Counter(concatenated_snps)
with open("./data/FISH_560/passing_HWE.snps", 'w') as OUTFILE:
    for snp, count in snp_counts.items():
        if count >= HWE_POP_THRESHOLD:
            OUTFILE.write("{}\n".format(snp))
        

print(len(concatenated_snps))  