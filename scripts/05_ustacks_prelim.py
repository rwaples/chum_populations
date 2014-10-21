import os
import os.path
import glob

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




fqgz_dir = "/media/Shared/Data/chum/populations/cleanSeqs"
ustacks_cmd = "ustacks -f {} -o /media/Shared/Data/chum/populations/stacks -i {} --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq "


for count, xx in enumerate(catalog_individuals):
    ind_id = count + 1
    ind_filename = os.path.join(fqgz_dir, xx + ".fq.gz")
    print(ustacks_cmd.format(ind_filename, ind_id))
    

other_individuals = ([os.path.basename(xx)[:-6] for xx in glob.glob(os.path.join(fqgz_dir, "*.fq.gz")) if os.path.basename(xx)[:-6] not in catalog_individuals])

for count, xx in enumerate(other_individuals):
    ind_id = count + 22
    ind_filename = os.path.join(fqgz_dir, xx + ".fq.gz")
    print(ustacks_cmd.format(ind_filename, ind_id))

ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0001.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 1 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 2 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 3 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 4 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 5 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0016.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 6 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0017.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 7 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0030.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 8 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0040.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 9 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0040.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 10 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 11 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0038.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 12 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 13 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0047.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 14 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0044.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 15 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0099.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 16 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0037.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 17 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0075.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 18 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0032.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 19 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0090.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 20 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 21 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 

ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 22 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 23 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 24 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 25 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0019.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 26 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0031.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 27 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 28 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0029.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 29 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0033.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 30 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0056.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 31 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0044.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 32 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0046.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 33 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0048.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 34 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0049.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 35 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0050.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 36 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0026.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 37 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0034.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 38 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0002.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 39 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0012.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 40 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0026.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 41 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0004.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 42 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0041.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 43 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 44 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0035.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 45 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0077.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 46 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0079.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 47 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0080.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 48 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0085.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 49 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0087.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 50 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0089.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 51 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0095.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 52 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0096.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 53 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0033.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 54 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0035.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 55 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0036.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 56 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0029.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 57 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0033.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 58 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0036.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 59 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0037.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 60 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0043.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 61 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0055.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 62 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0059.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 63 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0093.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 64 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0095.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 65 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0076.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 66 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0083.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 67 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0084.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 68 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0085.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 69 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0032.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 70 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0033.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 71 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0034.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 72 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0047.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 73 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0067.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 74 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0080.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 75 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0002.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 76 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0003.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 77 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0004.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 78 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0007.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 79 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 80 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 81 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0012.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 82 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0015.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 83 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMUW10_0016.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 84 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0014.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 85 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0015.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 86 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0016.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 87 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0017.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 88 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0018.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 89 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0022.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 90 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0024.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 91 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0025.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 92 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0005.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 93 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0006.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 94 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0007.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 95 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 96 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0010.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 97 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0012.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 98 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0017.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 99 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0024.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 100 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0025.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 101 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0036.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 102 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0037.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 103 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0030.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 104 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0031.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 105 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0032.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 106 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0047.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 107 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0048.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 108 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0053.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 109 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0061.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 110 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMLILLIW11_0063.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 111 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0012.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 112 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0018.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 113 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0019.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 114 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0022.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 115 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0024.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 116 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0026.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 117 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0027.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 118 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0005.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 119 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 120 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 121 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0045.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 122 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0046.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 123 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0047.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 124 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0052.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 125 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0058.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 126 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0059.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 127 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0061.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 128 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0064.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 129 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0066.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 130 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0074.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 131 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0091.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 132 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0001.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 133 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0003.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 134 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 135 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0016.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 136 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0017.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 137 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0018.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 138 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0019.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 139 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0020.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 140 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0022.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 141 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0053.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 142 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0056.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 143 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0057.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 144 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0059.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 145 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0028.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 146 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0039.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 147 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0056.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 148 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94F_0086.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 149 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0015.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 150 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0050.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 151 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0062.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 152 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0087.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 153 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0002.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 154 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 155 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0041.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 156 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0088.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 157 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0106.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 158 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0005.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 159 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0036.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 160 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0006.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 161 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0020.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 162 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0061.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 163 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0066.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 164 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0068.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 165 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0071.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 166 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0075.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 167 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0081.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 168 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0083.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 169 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSHERW94S_0085.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 170 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0004.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 171 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0006.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 172 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0007.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 173 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0008.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 174 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0011.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 175 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0012.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 176 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0013.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 177 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0018.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 178 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0033.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 179 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0035.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 180 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSKOO10_0036.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 181 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0058.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 182 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0062.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 183 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0069.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 184 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0072.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 185 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0085.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 186 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0089.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 187 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0090.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 188 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0102.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 189 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0104.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 190 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSNOH10_0105.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 191 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0001.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 192 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0004.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 193 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0038.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 194 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0046.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 195 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0052.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 196 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0055.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 197 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSQUA10_0058.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 198 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0009.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 199 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
ustacks -f /media/Shared/Data/chum/populations/cleanSeqs/CMSTILL10_0014.fq.gz -o /media/Shared/Data/chum/populations/stacks -i 200 --model_type bounded --bound_high 0.05 --alpha 0.1  -m 2 -M 4 -H -r --max_locus_stacks 4 -p 4 -t gzfastq 
