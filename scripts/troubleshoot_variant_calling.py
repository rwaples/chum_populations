import glob 
import os

def parse_name_fqgz(ind_filename):
    basename = os.path.basename(ind_filename)
    root, ext = os.path.splitext(basename)
    sample_name, fq = os.path.splitext(root)
    silli, ind_num = sample_name.split('_')
    return(basename, sample_name, silli, ind_num)
    
files_to_align = glob.glob('/media/Shared/Data/chum/populations/troubleshoot_variant_calling/clean_seqs/*.fq.gz')
# names are changed by adding a 'c'
my_ref_file = '/home/ipseg/Desktop/waples/chum_populations/data/ref/batch_42_CURATED.fasta.txt'

# this ref padding with N's at the end of each locus
# THIS SEEMS TO BE IMPORTANT!!!!
my_ref_file = '/media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_CURATED.fasta.txt'
samtools faidx '/media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_CURATED.fasta.txt'
bowtie-build '/media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_CURATED.fasta.txt' batch_42_bowtie1
bowtie2-build '/media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_CURATED.fasta.txt' batch_42_bowtie2


bwa_cmd = "bwa mem -t 6 -R '{}' {} {} | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bwa/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    read_group = "@RG\\tID:{}\\tSM:{}".format(sample_name, sample_name)
    print(bwa_cmd.format(read_group, my_ref_file, individual_fq_gz, sample_name))
    
# bowtie2
# NOTICE, trims the 3' end by nine bps to match length in full
bowtie2_cmd = "gunzip -c {} | bowtie2 -p 6 --trim3 9 --rg-id {} --rg SM:{} --rg PL:ILLUMINA /media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_bowtie2 - | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    print(bowtie2_cmd.format(individual_fq_gz, sample_name, sample_name, sample_name))

                                                              
# bowtie1
bowtie1_cmd = "gunzip -c {} | bowtie -p 6 --trim3 10 -v 3 -S --sam-RG ID:{} --sam-RG SM:{} --sam-RG PL:ILLUMINA /media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref/batch_42_bowtie1 - | samtools view -bhu - | samtools sort -m 2G -O bam -T temp_prefix -@ 2 - > /media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie1/{}.bam"

for individual_fq_gz in files_to_align:
    basename, sample_name, silli, ind_num = parse_name_fqgz(individual_fq_gz)
    print(bowtie1_cmd.format(individual_fq_gz, sample_name, sample_name, sample_name))
# for bowtie1 RG must be added to each seq line, possible to do with bamaddrg program

#find BAM files
BAM_files = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie1/*.bam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/*.bam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bwa/*.bam")

for xx in BAM_files:
    print("samtools view -h -o {} {}".format(xx.replace('.bam', '.sam'), xx))
    
SAM_files = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie1/*.sam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/*.sam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bwa/*.sam")

# filter out reads not that do not start their alignment at POS=1
for xx in SAM_files:
    with open(xx) as INFILE:
        with open(xx.replace(".sam", ".pos1.sam"), 'w') as OUTFILE:
            with open(xx.replace(".sam", ".discards.sam"), 'w') as DISCARDS:
                for line in INFILE:
                    if line.startswith("@"):
                        OUTFILE.write(line)
                    else:
                        if line.split("\t")[3] == '1':
                            OUTFILE.write(line)
                        else:
                            DISCARDS.write(line)

pos1_SAM_files = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie1/*.pos1.sam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/*.pos1.sam") + glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bwa/*.pos1.sam")

for xx in pos1_SAM_files:
    print("samtools view -Sb  {}  >  {}".format(xx, xx.replace('.sam', '.bam')))
    print("samtools index {}".format(xx.replace('.sam', '.bam')))

pos1_BAM_files_bowtie1 = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie1/*.pos1.bam") 
pos1_BAM_files_bowtie2 = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/*.pos1.bam")
pos1_BAM_files_bwa = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bwa/*.pos1.bam")

"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
--min-repeat-size 10 --binomial-obs-priors-off \
-b " + " -b ".join(pos1_BAM_files_bowtie1) + " --vcf  /media/Shared/Data/chum/populations/troubleshoot_variant_calling/test_bowtie1.raw.vcf"



BAM_files_bowtie2 = glob.glob("/media/Shared/Data/chum/populations/troubleshoot_variant_calling/aln/bowtie2/*.bam")

"freebayes -f /media/Shared/Data/chum/populations/troubleshoot_variant_calling/ref2/batch_42_CURATED.fasta.txt \
--binomial-obs-priors-off \
--bam " + " --bam ".join(BAM_files_bowtie2) + " --vcf  /media/Shared/Data/chum/populations/troubleshoot_variant_calling/test_bowtie2.raw.vcf" 


"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
--min-repeat-size 10 --binomial-obs-priors-off \
-b " + " -b ".join(pos1_BAM_files_bwa) + " --vcf  /media/Shared/Data/chum/populations/troubleshoot_variant_calling/test_bwa.raw.vcf"









# freebayes command
"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
--min-repeat-size 10 \
--binomial-obs-priors-off \
--populations /home/ipseg/Desktop/waples/chum_populations/data/fb.populations \
--bam-list /media/Shared/Data/chum/populations/aln/curated/bowtie2/bam.list \
--vcf /media/Shared/Data/chum/populations/troubleshoot_variant_calling/all.bt2.raw.vcf"
