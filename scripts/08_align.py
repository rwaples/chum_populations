import glob 
import os

# build index with bwa
runthis = 'bwa index ./chum_ref_batch_01.fa'
runthis = 'bwa index ./chum_ref_batch_02.fa'

def parse_name_fqgz(ind_filename):
    basename = os.path.basename(ind_filename)
    root, ext = os.path.splitext(basename)
    sample_name, fq = os.path.splitext(root)
    silli, ind_num = sample_name.split('_')
    return(basename, sample_name, silli, ind_num)


files_to_align = glob.glob('/media/Shared/Data/chum/populations/cleanSeqs/*.fq.gz')

my_ref_file = '/home/ipseg/Desktop/waples/chum_populations/data/ref/batch_42_CURATED.fasta.txt'

# could try to apply a large clipping penalty "-L 1000"
# should we try to trim the 3' end by nine units?
bwa_cmd = "bwa mem -t 6 -R '{}' {} {} | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/aln/curated/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    read_group = "@RG\\tID:{}\\tSM:{}".format(sample_name, sample_name)
    print(bwa_cmd.format(read_group, my_ref_file, individual_fq_gz, sample_name))
    
        
# bowtie1
bowtie1_cmd = "gunzip -c {} | bowtie -p 6 --trim3 9 -v 3 -S --sam-RG ID:{} --sam-RG SM:{} /media/Shared/Data/chum/populations/aln/batch_02/ref/bt1 - | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/aln/batch_02/bowtie1/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    print(bowtie1_cmd.format(individual_fq_gz, sample_name, sample_name, sample_name))
# for bowtie1 RG must be added to each seq line, possible to do with bamaddrg program 


# bowtie2  
bowtie2_cmd = "gunzip -c {} | bowtie2 -p 6 --trim3 9 --no-unal --rg-id {} --rg SM:{} --rg PL:ILLUMINA /home/ipseg/Desktop/waples/chum_populations/data/ref/batch_42_CURATED - | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/aln/curated/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    print(bowtie2_cmd.format(individual_fq_gz, sample_name, sample_name, sample_name))
                    