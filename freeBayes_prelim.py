import glob
import os.path

# prep for freebayes

# index the reference
samtools faidx ./chum_ref_batch_02.fa

#generate commands to convert sams to bams
# NO LONGER NEEDED
sam_files = glob.glob('/media/Shared/Data/chum/populations/aln/*.sam')
for SAM in sam_files:
    print("samtools view -bhu {sam_file} | samtools sort -m 4G -o {bam_file} -O bam -T temp_prefix -@ 2 -".format(sam_file = SAM, bam_file = SAM.replace(".sam", ".bam")))


# index each bam
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie2/*.bam')
for BAM in bam_files:
    print("samtools index {bam_file}".format(bam_file=BAM))


# make a CNVmap to represent duplicated loci, base ploidy will be set a 2, duplicated loci will be set at 4.
# CNV variation must be set in each sample:
# example BED line ([tab delim]):
# reference_sequence, start, end, sample name, copy number
mapped_stats = "/media/Shared/Data/chum/populations/map/mapped.stats.txt"
#my_BED_file = "/media/Shared/Data/chum/populations/fb/batch_01.bed"
my_BED_file = "/media/Shared/Data/chum/populations/fb/batch_02.bed"
my_target_file = "/media/Shared/Data/chum/populations/fb/batch_02_targets.bed"
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie1/*.bam')
my_sample_names = [os.path.splitext(os.path.basename(xx))[0] for xx in bam_files]

def generate_CNV_BED(map_stats_file, BED_file, sample_names, catID_col = 0, best_model_col = 4):
    with open(map_stats_file) as INFILE:
        with open (BED_file, 'w') as OUTFILE:
            skip_first = next(INFILE)
            for line in INFILE:
                line_split = line.strip().split()
                catID = line_split[catID_col]
                #print(catID)
                model_call = line_split[best_model_col]
                if (model_call != "AB"):
                    #print(line_split[catID_col], line_split[best_model_col])
                    for sample_name in sample_names:
                        OUTFILE.write("catID|{}\t-1\t-1\t{}\t4\n".format(catID, sample_name))

                    
generate_CNV_BED(mapped_stats, my_BED_file, my_sample_names)

def generate_targets_file(map_stats_file, targets_file, sample_names):
    with open(map_stats_file) as INFILE:
        with open (targets_file, 'w') as OUTFILE:
            skip_first = next(INFILE)
            for line in INFILE:
                line_split = line.strip().split()
                catID = line_split[0]
                for sample_name in sample_names:
                    OUTFILE.write("catID|{}\t0\t94\t{}\n".format(catID, sample_name))

generate_targets_file(map_stats_file = mapped_stats, targets_file = my_target_file, sample_names = my_sample_names)

def generate_basic_BED(fa_ref, BED_output):
    with open(fa_ref) as INFILE:
        with open(BED_output, 'w') as OUTFILE:
            for line in INFILE:
                header = line.strip()
                seq = next(INFILE).strip()
                OUTFILE.write("{}\t0\t{}\n".format(header, len(seq)))

generate_basic_BED(fa_ref = "/media/Shared/Data/chum/populations/aln/chum_ref_batch_02.fa", BED_output = "/media/Shared/Data/chum/populations/fb/batch_02_nosample.bed")

my_pops = [xx.split("_")[0] for xx in my_sample_names]

populations_file  = "/media/Shared/Data/chum/populations/fb/batch_01.populations"
with open(populations_file, 'w') as OUTFILE:
    for sn, pop in  zip(my_sample_names, my_pops):
        OUTFILE.write("{}\t{}\n".format(sn,pop))
        
# to examine coverage
# 
./bedtools genomecov -bga -ibam '//media/Shared/Data/chum/populations/aln/batch_02/bowtie1/CMUW10X_0008.bam' > '/home/ipseg/Desktop/waples/chum_populations/basic_stats' 

# also we may want to increase the clipping penalty in bwa 

-D --read-dependence-factor N
                   Incorporate non-independence of reads by scaling successive
                   observations by this factor during data likelihood
                   calculations.  default: 0.9

-K --pooled-continuous
            Output all alleles which pass input filters, regardles of
               genotyping outcome or model.


# vcf flags to check?
##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality">

piping to bowtie, noticed that old chum reads were trimmed to 85 bases 
# gunzip -c /media/Shared/Data/chum/populations/cleanSeqs/CMUW10X_0008.fq.gz | bowtie -3 9 -S --sam-nohead -k 2 bt1 -  > '/media/Shared/Data/chum/populations/work/chum_08.sam'

# construct command to run freebayes
parental_bam_files = ["/media/Shared/Data/chum/populations/aln/batch_02/CMUW10X_0001.bam", "/media/Shared/Data/chum/populations/aln/CMUW10X_0008.bam", "/media/Shared/Data/chum/populations/aln/CMUW10X_0009.bam"]
test_bam_files = ["/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0030.bam", "/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0033.bam", "/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0040.bam"]

bam_files = glob.glob('/media/Shared/Data/chum/populations/aln//batch_02/*X*.bam')
bam_files += glob.glob('/media/Shared/Data/chum/populations/aln//batch_02/CMSHER*.bam')

bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie1/new.bam')
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie2/*.bam')


"freebayes -f /media/Shared/Data/chum/populations/aln/batch_02/ref/chum_ref_batch_02.fa \
--min-coverage 2 \
--min-repeat-size 100 \
--haplotype-length 0 \
--min-alternate-fraction 0.05 \
--min-alternate-count 1 \
--min-alternate-total 2 \
--binomial-obs-priors-off \
--theta .01 \
--targets /media/Shared/Data/chum/populations/fb/batch_02_nosample.bed \
--cnv-map /media/Shared/Data/chum/populations/fb/batch_02.bed \
-m 10 \
-q 10 \
-b " + " -b ".join(bam_files) + " --vcf /home/ipseg/Desktop/waples/chum_populations/batch_02_bt2.vcf"




--pooled-continuous \

--populations /media/Shared/Data/chum/populations/fb/batch_01.populations \


"freebayes -f /media/Shared/Data/chum/populations/aln/batch_02/chum_ref_batch_02.fa \
--min-alternate-fraction 0.1 \
-e 2 \
--ploidy 4 \
--exclude-unobserved-genotypes \
--genotype-qualities \
-b " + " -b ".join(bam_files) + " --vcf ./batch_03.vcf"


-q 10 \
-m 10 \





TODO

use --populations to separate my populations
--min-alternate-fraction "should be changed for ploidy > 2"

examine 

--no-population-priors 
--pooled-discrete
--hwe-priors-off
--cnv-map
	BED file =  reference sequence, start, end, sample name, copy number
