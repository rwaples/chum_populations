# pstacks
import os
import os.path
import glob

bam_dir = "/media/Shared/Data/chum/populations/aln/curated/bowtie2/start_filter"

bam_dir = "/media/Shared/Data/chum/populations/aln/curated/bwa"

pstacks_cmd = "pstacks -f {} -o /media/Shared/Data/chum/populations/pstacks -i {} --model_type bounded --bound_high 0.05 --alpha 0.05 -m 2 -p 6 -t bam"

sstacks_cmd = "sstacks -b 2 -c batch_2 -o /media/Shared/Data/chum/populations/pstacks -p 6 -g -s {} "

#ind_bams = ([os.path.basename(xx)[:-4] for xx in glob.glob(os.path.join(bam_dir, "*.bam"))])

ind_bams = glob.glob(os.path.join(bam_dir, "*.bam"))


for count, xx in enumerate(ind_bams):
    ind_id = count + 1
    ind_filename = os.path.join(bam_dir, "{}".format(xx))
    print(pstacks_cmd.format(ind_filename, ind_id))


for xx in ind_bams:
    yy = os.path.basename(xx)[:-4]
    print(sstacks_cmd.format(yy))



