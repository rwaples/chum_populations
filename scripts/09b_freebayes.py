import os
import glob


bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/curated/bowtie2/*.bam')


# freebayes command
"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED_with_Ns.fasta \
--min-mapping-quality 10 \
--min-base-quality 10 \
--min-repeat-size 10 \
--binomial-obs-priors-off \
--bam /media/Shared/Data/chum/populations/aln/curated/bowtie2/CMSTILL10_0014.bam \
--vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.raw.vcf"


#--bam-list /media/Shared/Data/chum/populations/aln/curated/bowtie2/bam.list \



--populations /home/ipseg/Desktop/waples/chum_populations/data/fb.populations \
--cnv-map /home/ipseg/Desktop/waples/chum_populations/data/ref/bed/CNV.bed \
--bam-list /media/Shared/Data/chum/populations/aln/curated/bowtie2/bam.list \
--samples /home/ipseg/Desktop/waples/chum_populations/data/fb.all.samples \




# targets file currently has issues, doen't match with ref

--targets /home/ipseg/Desktop/waples/chum_populations/data/ref/bed/mapped.bed \




-D --read-dependence-factor N
                   Incorporate non-independence of reads by scaling successive
                   observations by this factor during data likelihood
                   calculations.  default: 0.9

-K --pooled-continuous
            Output all alleles which pass input filters, regardles of
               genotyping outcome or model.

piping to bowtie, noticed that old chum reads were trimmed to 85 bases 
# gunzip -c /media/Shared/Data/chum/populations/cleanSeqs/CMUW10X_0008.fq.gz | bowtie -3 9 -S --sam-nohead -k 2 bt1 -  > '/media/Shared/Data/chum/populations/work/chum_08.sam'

# construct command to run freebayes
parental_bam_files = ["/media/Shared/Data/chum/populations/aln/batch_02/CMUW10X_0001.bam", "/media/Shared/Data/chum/populations/aln/CMUW10X_0008.bam", "/media/Shared/Data/chum/populations/aln/CMUW10X_0009.bam"]
test_bam_files = ["/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0030.bam", "/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0033.bam", "/media/Shared/Data/chum/populations/aln/batch_02/CMHAMM10_0040.bam"]

bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/curated/bowtie2/*.bam')
#bam_files += glob.glob('/media/Shared/Data/chum/populations/aln//batch_02/CMSHER*.bam')



# all pops
# look at /home/ipseg/Programs/freebayes/scripts/samples.cnv for example cnv variation
#first without phasing 
"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-m 10 -q 10 --binomial-obs-priors-off \
--populations /media/Shared/Data/chum/populations/fb/populations \
-b " + " -b ".join(bam_files) + " --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/all.raw.vcf"



# freebayes-parrallel
"""/home/ipseg/Programs/freebayes/scripts/fasta_generate_regions_RW.py /home/ipseg/Desktop/waples/chum_populations/data/ref/batch_42_CURATED.fasta.txt 1000 > /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/regions.txt"""

# regions file is not cooperating - naming of chromosomes is not occuring correctly
"""/home/ipseg/Programs/freebayes/scripts/freebayes-parallel /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/freebayes/regions.txt 4 \
-f /home/ipseg/Desktop/waples/chum_populations/data/ref/batch_42_CURATED.fasta.txt /media/Shared/Data/chum/populations/aln/curated/CMUW10_0008.bam"""


bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie1/new.bam')
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/batch_02/bowtie2/*.bam')


"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
--genotype-qualities \
--binomial-obs-priors-off \
--cnv-map /media/Shared/Data/chum/populations/fb/batch_42_CURATED.bed \
-b " + " -b ".join(bam_files) + " --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_curated.vcf"

"freebayes -f /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
--genotype-qualities \
--min-coverage 2 \
--min-repeat-size 20 \
--haplotype-length 0 \
--min-alternate-fraction 0.05 \
--min-alternate-count 1 \
--min-alternate-total 2 \
--binomial-obs-priors-off \
--theta .01 \
--populations /media/Shared/Data/chum/populations/fb/populations \
--cnv-map /media/Shared/Data/chum/populations/fb/batch_42_CURATED.bed \
-D .85 \
-m 20 \
-q 10 \
-R 20 \
-U 5 \
-e 1 \
-b " + " -b ".join(bam_files) + " --vcf /home/ipseg/Desktop/waples/chum_populations/results/batch_42_curated.vcf"


"freebayes -f /media/Shared/Data/chum/populations/aln/batch_02/chum_ref_batch_02.fa \
--min-alternate-fraction 0.1 \
-e 2 \
--ploidy 4 \
--exclude-unobserved-genotypes \
--genotype-qualities \
-b " + " -b ".join(bam_files) + " --vcf ./batch_03.vcf"


--targets /media/Shared/Data/chum/populations/fb/batch_42_CURATED_targets.bed \ 