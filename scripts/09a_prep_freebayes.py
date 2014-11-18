import glob
import os.path

# prep for freebayes

# index the reference
samtools faidx /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt

#generate commands to convert sams to bams
# NO LONGER NEEDED
#sam_files = glob.glob('/media/Shared/Data/chum/populations/aln/*.sam')
#for SAM in sam_files:
#    print("samtools view -bhu {sam_file} | samtools sort -m 4G -o {bam_file} -O bam -T temp_prefix -@ 2 -".format(sam_file = SAM, bam_file = SAM.replace(".sam", ".bam")))

# index each bam
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/curated/*.bam')
for BAM in bam_files:
    print("samtools index {bam_file}".format(bam_file=BAM))


# make a CNVmap to represent duplicated loci, base ploidy will be set a 2, duplicated loci will be set at 4.
# CNV variation must be set in each sample:
# example BED line ([tab delim]):
# reference_sequence, start, end, sample name, copy number
mapped_stats = "/media/Shared/Data/chum/populations/map/mapped.stats.txt"
#my_BED_file = "/media/Shared/Data/chum/populations/fb/batch_01.bed"
my_BED_file = "/media/Shared/Data/chum/populations/fb/batch_42_CURATED.bed"
my_target_file = "/media/Shared/Data/chum/populations/fb/batch_42_CURATED_targets.bed"
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/curated/*.bam')
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
                        OUTFILE.write("{}\t-1\t-1\t{}\t4\n".format(catID, sample_name))

                    
generate_CNV_BED(mapped_stats, my_BED_file, my_sample_names)

def generate_targets_file(map_stats_file, targets_file, sample_names):
    with open(map_stats_file) as INFILE:
        with open (targets_file, 'w') as OUTFILE:
            skip_first = next(INFILE)
            for line in INFILE:
                line_split = line.strip().split()
                catID = line_split[0]
                for sample_name in sample_names:
                    OUTFILE.write("{}\t0\t84\t{}\n".format(catID, sample_name))

generate_targets_file(map_stats_file = mapped_stats, targets_file = my_target_file, sample_names = my_sample_names)

def generate_basic_BED(fa_ref, BED_output):
    with open(fa_ref) as INFILE:
        with open(BED_output, 'w') as OUTFILE:
            for line in INFILE:
                header = line.strip()
                seq = next(INFILE).strip()
                OUTFILE.write("{}\t0\t{}\n".format(header, len(seq)))

generate_basic_BED(fa_ref = "/media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt", BED_output = "/media/Shared/Data/chum/populations/fb/batch_42_CURATED_nosample.bed")

my_pops = [xx.split("_")[0] for xx in my_sample_names]

populations_file  = "/media/Shared/Data/chum/populations/fb/populations"
with open(populations_file, 'w') as OUTFILE:
    for sn, pop in  zip(my_sample_names, my_pops):
        OUTFILE.write("{}\t{}\n".format(sn,pop))
        








TODO

use --populations to separate my populations
--min-alternate-fraction "should be changed for ploidy > 2"
