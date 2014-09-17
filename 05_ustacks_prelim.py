
catalog_individuals = [
    # mappings parents
    'CMUW10X_0001', 'CMUW10X_0008', 'CMUW10X_0009',
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

fqgz_dir = "/media/Shared/Data/chum/populations/cleanSeqs"

for count, xx in enumerate(catalog_individuals):
    ind_id = count + 1
    ind_filename = os.path.join(fqgz_dir, xx + ".fq.gz")
    ustacks_cmd = "ustacks -f {} -o /media/Shared/Data/chum/populations/stacks -i {} -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq"
    print(ustacks_cmd.format(ind_filename, ind_id))
    

ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10X_0001.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 1 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10X_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 2 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10X_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 3 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 4 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 5 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0016.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 6 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0017.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 7 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0030.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 8 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0040.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 9 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0040.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 10 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 11 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0038.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 12 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 13 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0047.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 14 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0044.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 15 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0099.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 16 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0037.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 17 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0075.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 18 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0032.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 19 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0090.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 20 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 21 -m 3 -M 4 -H -r --max_locus_stacks 4 -p 6 -t gzfastq

