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
    '/media/Shared/Data/chum/populations/barcodes/CMX3_individuals.txt' : '/media/Shared/Data/chum/populations/cleanSeqs/CMX3'
    }

for rename_file, fq_dir in rename_dict.items():
    with open(rename_file) as INFILE:
        for line in INFILE:
            barcode, ind_name = line.strip().split('\t')
            #print os.path.exists((os.path.join(fq_dir, 'sample_{}.fq'.format(barcode))))
            old_name = os.path.join(fq_dir, 'sample_{}.fq'.format(barcode))
            new_name = os.path.join(fq_dir, '{}.fq'.format(ind_name))
            print(old_name, new_name)          
            os.rename(old_name, new_name)
            with open(os.path.join(fq_dir, 'rename_record.txt'), 'a') as OUTFILE:
                OUTFILE.write(old_name)
                OUTFILE.write("\t")
                OUTFILE.write(new_name)
                OUTFILE.write("\n")
            #with open(os.path.join(fq_dir, 'sample_{}.fq'.format(barcode))) as FQ_FILE:
            #   with gzip.open(os.path.join(fa_dir, '{}.fq.gz'.format(ind_name)), 'wb') as GZ_FILE:
            #       GZ_FILE.writelines(FQ_FILE)

PE_rename_dict = { 
    # path_to_renaming_file : path_to_target_dir
    '/media/Shared/Data/chum/PE/barcodes/chum_PE_individuals.txt' : '/media/Shared/Data/chum/PE/clean_demultiplex'}
    
for rename_file, fq_dir in PE_rename_dict.items():
    with open(rename_file) as INFILE:
        for line in INFILE:
            barcode, ind_name = line.strip().split('\t')
            #print os.path.exists((os.path.join(fq_dir, 'sample_{}.fq'.format(barcode))))
            old_name_1 = os.path.join(fq_dir, 'sample_{}.1.fq'.format(barcode))
            old_name_2 = os.path.join(fq_dir, 'sample_{}.2.fq'.format(barcode))
            new_name_1 = os.path.join(fq_dir, '{}.1.fq'.format(ind_name))
            new_name_2 = os.path.join(fq_dir, '{}.2.fq'.format(ind_name))
            print(old_name_1, old_name_2, new_name_1, new_name_2)          
            os.rename(old_name_1, new_name_1)
            os.rename(old_name_2, new_name_2)
            with open(os.path.join(fq_dir, 'rename_record.txt'), 'a') as OUTFILE:
                OUTFILE.write("\t".join([old_name_1, old_name_2, new_name_1, new_name_2]))
                OUTFILE.write("\n")
            #with open(os.path.join(fq_dir, 'sample_{}.fq'.format(barcode))) as FQ_FILE:
            #   with gzip.open(os.path.join(fa_dir, '{}.fq.gz'.format(ind_name)), 'wb') as GZ_FILE:
            #       GZ_FILE.writelines(FQ_FILE)            
            
