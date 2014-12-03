import glob
import os.path

# index the reference
shell command: samtools faidx /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt

# index each bam
shell command:  for f in *.bam ; do samtools index "$f"; done

# create a file listing each bam
shell command: ls -d $PWD/*.bam > ./bam.list


def generate_CNV_BED(map_stats_file, BED_file, sample_names, catID_col = 0, best_model_col = 4):
    """
    this file is passed to freebayes to specify tetrasomically-inherited loci
    format:
    reference [\t] start [\t] end [\t] sample name [\t] copy number
    "-1" for start + end specifies entire reference
    """
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

mapped_stats = "/home/ipseg/Desktop/waples/chum_populations/data/mapped.stats.txt"
my_CNV_BED_file = "/home/ipseg/Desktop/waples/chum_populations/data/ref/bed/CNV.bed"
bam_files = glob.glob('/media/Shared/Data/chum/populations/aln/curated/bowtie2/*.bam')
my_sample_names = [os.path.splitext(os.path.basename(xx))[0] for xx in bam_files]
generate_CNV_BED(mapped_stats, my_CNV_BED_file, my_sample_names)

def generate_mapped_targets_file(map_stats_file, targets_file):
    with open(map_stats_file) as INFILE:
        with open (targets_file, 'w') as OUTFILE:
            skip_first = next(INFILE)
            for line in INFILE:
                line_split = line.strip().split()
                catID = line_split[0]
                OUTFILE.write("{}\t0\t84\n".format(catID))

my_target_file = "/home/ipseg/Desktop/waples/chum_populations/data/ref/bed/mapped.bed"
generate_mapped_targets_file(map_stats_file = mapped_stats, targets_file = my_target_file)

my_pops = [xx.split("_")[0] for xx in my_sample_names]
populations_file  = "/home/ipseg/Desktop/waples/chum_populations/data/fb.populations"
with open(populations_file, 'w') as OUTFILE:
    for sn, pop in  zip(my_sample_names, my_pops):
        OUTFILE.write("{}\t{}\n".format(sn,pop))



all_samples_file  = "/home/ipseg/Desktop/waples/chum_populations/data/fb.all.samples"
CMUW_samples_file  = "/home/ipseg/Desktop/waples/chum_populations/data/fb.CMUW.samples"
with open(all_samples_file, 'w') as OUTFILE:
    for sn in  my_sample_names:
        OUTFILE.write("{}\n".format(sn))

with open(CMUW_samples_file, 'w') as OUTFILE:
    for sn in my_sample_names:
        OUTFILE.write("{}\n".format(sn))

#####

#generate commands to convert sams to bams
# NO LONGER NEEDED - done during alignment step
#sam_files = glob.glob('/media/Shared/Data/chum/populations/aln/*.sam')
#for SAM in sam_files:
#    print("samtools view -bhu {sam_file} | samtools sort -m 4G -o {bam_file} -O bam -T temp_prefix -@ 2 -".format(sam_file = SAM, bam_file = SAM.replace(".sam", ".bam")))


def generate_basic_BED(fa_ref, BED_output):
    with open(fa_ref) as INFILE:
        with open(BED_output, 'w') as OUTFILE:
            for line in INFILE:
                header = line.strip()
                seq = next(INFILE).strip()
                OUTFILE.write("{}\t0\t{}\n".format(header, len(seq)))

generate_basic_BED(fa_ref = "/media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt", BED_output = "/media/Shared/Data/chum/populations/fb/batch_42_CURATED_nosample.bed")



        
