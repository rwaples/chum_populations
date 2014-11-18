import pandas as pd
import matplotlib.pyplot as plt

# stats files generated with something like:
# bedtools genomecov -bga -ibam '/media/Shared/Data/chum/populations/aln/curated/CMSHERW94F_0010.bam' > '/home/ipseg/Desktop/waples/chum_populations/basic_stats' 


bowtie2_file = "/home/ipseg/Desktop/waples/chum_populations/basic_stats"
bwa_stats_file = "/home/ipseg/Desktop/waples/chum_populations/basic_stats2"

bt2 = pd.read_table(bowtie2_file, header = None, index_col = 0)
bwa = pd.read_table(bwa_stats_file, header = None, index_col = 0)

bt2_max = bt2.groupby(level =0).max()
bwa_max = bwa.groupby(level =0).max()
bt2_max.columns = ['srt', 'stop', 'depth']
bwa_max.columns = ['srt', 'stop', 'depth']

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
