import glob
import os
import collections
import sys

output_file = sys.argv[1]
snplist_files = sys.argv[2:]

# build a list of all the snps from each family.  
# snps will be present x times if they appear in x families
concatenated_snps = list()
for snplist_file in snplist_files:
    with open(snplist_file, 'r') as INFILE:
        snps = [xx.strip() for xx in INFILE.readlines()]
        concatenated_snps += snps
        
# write any snps that appear at least once
with open(output_file, 'w') as OUTFILE:
    for snp in set(concatenated_snps):
        OUTFILE.write("{}\n".format(snp))
        
print(len(set(concatenated_snps)))