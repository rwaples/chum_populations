catalog_individuals = [
# mappings parents
'CMUW10_0001', 'CMUW10_0008', 'CMUW10_0009',
# CMKALA03
'CMKALA03_0009', 'CMKALA03_0011',
# CMSQUA10
'CMSQUA10_0016', 'CMSQUA10_0017',
# CMHAMM10
'CMHAMM10_0030', 'CMHAMM10_0040',
# CMSKOO10
'CMSKOO10_0040', 'CMSKOO10_0010',
# CMSNOH10
'CMSNOH10_0038', 'CMSNOH10_0027',
# CMSTILL10
'CMSTILL10_0047', 'CMSTILL10_0044',
# CMLILLIW11
'CMLILLIW11_0099', 'CMLILLIW11_0037',
# CMSHERW94F
'CMSHERW94F_0075', 'CMSHERW94F_0032',
# CMSHERW94S
'CMSHERW94S_0090', 'CMSHERW94S_0027'
]

print('cstacks -b 3  -o /media/Shared/Data/chum/populations/pstacks/batch_3 -p 6 -g' + " -s ./" + " -s ./".join(catalog_individuals))

#remove parents from catalog indiviudals and instead use the pre-existing catalog from the mappping paper, after convserion to new stacks format
# use /home/ipseg/Programs/Stacks/stacks-1.20/scripts/convert_stacks.pl
# from stacks directory, run this
#print('cstacks -b 03 --catalog /media/Shared/Data/chum/populations/stacks/mapping_catalog/batch_42' + " -s ./" + " -s ./".join(catalog_individuals) +' -o ./ -n 4 -p 6')

print("cstacks -b 05 -s ./" + " -s ./".join(catalog_individuals) + " -o ./ -n 4 -p 6")


# sstacks _prelim
for count, xx in enumerate(catalog_individuals):
    ind_id = count + 1
    sstacks_cmd = "sstacks -b 10 -c ./batch_10 -s ./{} -o /media/Shared/Data/chum/populations/stacks -p 8"
    print(sstacks_cmd.format(xx))
    
    
# sstacks _prelim
for count, xx in enumerate(other_individuals):
    ind_id = count + 1
    sstacks_cmd = "sstacks -b 10 -c ./batch_10 -s ./{} -o /media/Shared/Data/chum/populations/stacks -p 8"
    print(sstacks_cmd.format(xx))
    
    
    
    
genotypes -b 10 -P /media/Shared/Data/chum/populations/stacks -s 
populations -b 10 -P /media/Shared/Data/chum/populations/stacks -s -t 6 -r .5 -p .5 -a .05 --fstats --genepop --plink 
populations -b 10 -P /media/Shared/Data/chum/populations/stacks -s -t 6 -r .5 -p 4 -a .05 --fastphase --phase --beagle_phase --genepop --plink --phylip_var --vcf -M /media/Shared/Data/chum/populations/stacks/pop_map/pop_map.txt -W /media/Shared/Data/chum/populations/stacks_output/populations/whitelist.txt

# rxstacks_prelim
rxstacks -b 10 -P ./ -o ./rxstacks_b10 --model_type bounded --bound_high 0.05 --lnl_filter --lnl_dist --lnl_lim -10.0 -t 6
rxstacks -b 3 -P ./ -o ./rxstacks_2 --model_type bounded --bound_high 0.05 --lnl_filter --lnl_lim -10.0 --lnl_dist --max_haplo 3 --prune_haplo -t 6



sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0010 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0008 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0009 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0010 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0019 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0031 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0027 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0029 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0033 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0056 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0044 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0046 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0048 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0049 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0050 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0026 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0034 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0002 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0012 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0026 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0004 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0041 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0011 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0035 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0077 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0079 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0080 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0085 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0087 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0089 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0095 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0096 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0033 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0035 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0036 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0029 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0033 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0036 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0037 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0043 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0055 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0059 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0093 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0095 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0076 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0083 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0084 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0085 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0032 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0033 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0034 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0047 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0067 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0080 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0002 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0003 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0004 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0007 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0010 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0011 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0012 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0015 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMUW10_0016 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0014 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0015 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0016 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0017 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0018 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0022 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0024 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0025 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0005 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0006 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0007 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0008 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0010 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0012 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0017 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0024 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0025 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0036 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMKALA03_0037 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0030 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0031 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0032 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0047 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0048 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0053 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0061 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMLILLIW11_0063 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0012 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0018 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0019 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0022 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0024 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0026 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0027 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0005 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0008 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMHAMM10_0011 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0045 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0046 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0047 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0052 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0058 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0059 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0061 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0064 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0066 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0074 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0091 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0001 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0003 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0009 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0016 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0017 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0018 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0019 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0020 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0022 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0053 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0056 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0057 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0059 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0028 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0039 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0056 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94F_0086 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0015 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0050 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0062 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0087 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0002 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0009 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0041 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0088 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0106 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0005 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0036 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0006 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0020 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0061 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0066 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0068 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0071 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0075 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0081 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0083 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSHERW94S_0085 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0004 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0006 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0007 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0008 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0011 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0012 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0013 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0018 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0033 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0035 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSKOO10_0036 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0058 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0062 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0069 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0072 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0085 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0089 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0090 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0102 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0104 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSNOH10_0105 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0001 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0004 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0038 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0046 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0052 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0055 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSQUA10_0058 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0009 -o /media/Shared/Data/chum/populations/stacks -p 8
sstacks -b 10 -c ./batch_10 -s ./CMSTILL10_0014 -o /media/Shared/Data/chum/populations/stacks -p 8
