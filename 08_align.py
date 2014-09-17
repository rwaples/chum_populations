import glob 
import os

# build index with bwa
runthis = 'bwa index ./chum_ref_batch_01.fa'

files_to_align = glob.glob('/media/Shared/Data/chum/populations/cleanSeqs/*.fq.gz')
my_ref_file = '/media/Shared/Data/chum/populations/aln/chum_ref_batch_01.fa'


def parse_name_fqgz(ind_filename):
    basename = os.path.basename(ind_filename)
    root, ext = os.path.splitext(basename)
    sample_name, fq = os.path.splitext(root)
    silli, ind_num = sample_name.split('_')
    return(basename, sample_name, silli, ind_num)

bwa_cmd = "bwa mem -t 6 -R '{}' {} {} > /media/Shared/Data/chum/populations/aln/{}.sam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    read_group = "@RG\\tID:{}\\tSM:{}".format(sample_name, sample_name)
    print(bwa_cmd.format(read_group, my_ref_file,  individual_fq_gz, sample_name))
    
        

                    