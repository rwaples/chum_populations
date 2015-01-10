import pandas as pd
import matplotlib.pyplot as plt


# stats files generated with something like:
# bedtools genomecov -bga -ibam '/media/Shared/Data/chum/populations/aln/curated/CMSHERW94F_0010.bam' > '/home/ipseg/Desktop/waples/chum_populations/basic_stats' 

# Use bedtools to summarize coverage across each locus
"""
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bowtie2/CMSHERW94S_0001.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMSHERW94S_0001.bowtie2.txt
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bowtie2/CMLILLIW11_0032.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0032.bowtie2.txt
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bowtie2/CMLILLIW11_0079.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0079.bowtie2.txt
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bwa/CMSHERW94S_0001.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMSHERW94S_0001.bwa.txt
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bwa/CMLILLIW11_0032.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0032.bwa.txt
bedtools genomecov -bga -ibam /media/Shared/Data/chum/populations/aln/curated/bwa/CMLILLIW11_0079.bam > /home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0079.bwa.txt 
"""

CMSHERW94S_0001_bowtie2_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMSHERW94S_0001.bowtie2.txt"
CMLILLIW11_0032_bowtie2_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0032.bowtie2.txt"
CMLILLIW11_0079_bowtie2_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0079.bowtie2.txt"
CMSHERW94S_0001_bwa_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMSHERW94S_0001.bwa.txt"
CMLILLIW11_0032_bwa_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0032.bwa.txt"
CMLILLIW11_0079_bwa_filename = "/home/ipseg/Desktop/waples/chum_populations/data/aln_stats/CMLILLIW11_0079.bwa.txt"

CMSHERW94S_0001_bt2 = pd.read_table(CMSHERW94S_0001_bowtie2_filename, header = None, index_col = 0)
CMLILLIW11_0032_bt2 = pd.read_table(CMLILLIW11_0032_bowtie2_filename, header = None, index_col = 0)
CMLILLIW11_0079_bt2 = pd.read_table(CMLILLIW11_0079_bowtie2_filename, header = None, index_col = 0)
CMSHERW94S_0001_bwa = pd.read_table(CMSHERW94S_0001_bwa_filename, header = None, index_col = 0)
CMLILLIW11_0032_bwa = pd.read_table(CMLILLIW11_0032_bwa_filename, header = None, index_col = 0)
CMLILLIW11_0079_bwa = pd.read_table(CMLILLIW11_0079_bwa_filename, header = None, index_col = 0)



CMSHERW94S_0001_bt2.columns = ['srt', 'stop', 'depth']
CMLILLIW11_0032_bt2.columns = ['srt', 'stop', 'depth']
maxes = pd.merge(CMSHERW94S_0001_bt2, CMLILLIW11_0032_bt2, left_index = True, right_index = True, how = 'outer')
plt.plot(maxes.depth_x, maxes.depth_y, 'bo', alpha = .2)
pyplot.yscale('log')
pyplot.xscale('log')
pyplot.ylim(.8, 12000)
pyplot.xlim(.8, 12000)



my_data_tables = [CMSHERW94S_0001_bt2, CMLILLIW11_0032_bt2, CMLILLIW11_0079_bt2, CMSHERW94S_0001_bwa, CMLILLIW11_0032_bwa, CMLILLIW11_0079_bwa]

def find_depth(data_table):
    maxes_per_locus = data_table.groupby(level =0).max()
    maxes_per_locus.columns = ['srt', 'stop', 'depth']
    return(maxes_per_locus)

my_depths = [find_depth(xx) for xx in my_data_tables]


merge_1 = pd.merge(my_depths[0], my_depths[1], left_index = True, right_index = True, how = 'outer')
merge_2 = pd.merge(merge_1, my_depths[2], left_index = True, right_index = True, how = 'outer')
merge_3 = pd.merge(merge_2, my_depths[3], left_index = True, right_index = True, how = 'outer')
merge_4 = pd.merge(merge_3, my_depths[4], left_index = True, right_index = True, how = 'outer')
depths = pd.merge(merge_4, my_depths[5], left_index = True, right_index = True, how = 'outer')

pd.concat(my_depths, join = "outer", verify_integrity = True)

bt2_sum = np.sum(bt2_max.depth)
bwa_sum = np.sum(bwa_max.depth)

bt2_max.depth[bt2_max.depth >= 10000] = 10000
bwa_max.depth[bwa_max.depth >= 10000] = 10000

maxes = pd.merge(bt2_max, bwa_max, left_index = True, right_index = True, how = 'outer')
plt.plot(maxes.depth_x, maxes.depth_y, 'bo', alpha = .2)
pyplot.yscale('log')
pyplot.ylim(.8, 12000)
pyplot.xscale('log')
pyplot.xlim(.8, 12000)



# OLD
#bowtie2_file = "/home/ipseg/Desktop/waples/chum_populations/basic_stats"
#bwa_stats_file = "/home/ipseg/Desktop/waples/chum_populations/basic_stats2"