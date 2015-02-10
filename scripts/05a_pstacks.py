# pstacks
import os
import os.path
import glob


bam_dir = "/media/Shared/Data/chum/populations/aln/curated/bowtie2"



pstacks_cmd = "pstacks -f {} -o /media/Shared/Data/chum/populations/pstacks -i {} --model_type bounded --bound_high 0.05 --alpha 0.1  -m 1 -p 6 -t bam"

sstacks_cmd = "sstacks -b 1 -c batch_1 -o /media/Shared/Data/chum/populations/pstacks -s {} -p 6  -g"

ind_bams = ([os.path.basename(xx)[:-6] for xx in glob.glob(os.path.join(bam_dir, "*.bam"))])

for count, xx in enumerate(individuals):
    ind_id = count + 1
    ind_filename = os.path.join(bam_dir, "{}.bam".format(xx))
    print(pstacks_cmd.format(ind_filename, ind_id))



for xx in individuals:
    print(sstacks_cmd.format(xx))



