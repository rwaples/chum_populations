import os
import os.path
import gzip


#import glob
#glob.glob("/media/Shared/Data/chum/populations/barcodes/*_individuals.txt")


rename_dict = { 
    # path_to_renaming_file : path_to_target_dir
    '/media/Shared/Data/chum/populations/barcodes/CM01_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CM01',
    '/media/Shared/Data/chum/populations/barcodes/CM05_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CM05',
    '/media/Shared/Data/chum/populations/barcodes/CM06_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CM06',
    '/media/Shared/Data/chum/populations/barcodes/CM09_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CM09',
    '/media/Shared/Data/chum/populations/barcodes/CMX1_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CMX1',
#    '/media/Shared/Data/chum/populations/barcodes/CMX3_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CMX3'
    }



for rename_file, fa_dir in rename_dict.items():
    with open(rename_file) as INFILE:
        for line in INFILE:
            barcode, ind_name = line.strip().split('\t')
            #print os.path.exists((os.path.join(fa_dir, 'sample_{}.fa'.format(barcode))))           
            #os.rename(os.path.join(fa_dir, 'sample_{}.fa'.format(barcode)), os.path.join(fa_dir, '{}.fa'.format(ind_name)))
            with open(os.path.join(fa_dir, 'sample_{}.fa'.format(barcode))) as FA_FILE:
                with gzip.open(os.path.join(fa_dir, '{}.fa.gz'.format(ind_name)), 'wb') as GZ_FILE:
                    for line in FA_FILE:
                        GZ_FILE.write(line)


            
            
