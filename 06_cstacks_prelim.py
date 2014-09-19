catalog_individuals = [
    # mappings parents
    'CMUW10X_0001', 'CMUW10X_0008', 'CMUW10X_0009'],
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


#remove paretns from catalog indiviudals and instead use the pre-existing catalog from the mappping paper, after convserion to new stacks format
# use /home/ipseg/Programs/Stacks/stacks-1.20/scripts/convert_stacks.pl

# from stacks directory, run this
print('cstacks -b 03 --catalog /media/Shared/Data/chum/populations/stacks/mapping_catalog/batch_42' + " -s ./" + " -s ./".join(catalog_individuals) +' -o ./ -n -p 6')

print('cstacks -b 03 ' + " -s ./" + " -s ./".join(catalog_individuals) +' -o ./ -n -p 6')


#cstacks -b 01 -s **** -o /media/Shared/Data/chum/populations/stacks -n -p 6


