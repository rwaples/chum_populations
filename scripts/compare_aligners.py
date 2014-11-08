import pandas as pd
import matplotlib.pyplot as plt

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

bt2_max.depth[bt2_means.depth >= 1000] = 1000
bwa_max.depth[bwa_means.depth >= 1000] = 1000

means = pd.merge(bt2_means, bwa_means, left_index = True, right_index = True, how = 'outer')
plt.plot(means.depth_x, means.depth_y, 'bo')
pyplot.yscale('log')
pyplot.ylim(.8, 2000)
pyplot.xscale('log')
pyplot.xlim(.8, 2000)

measn.ix('3_x')
set(bt2.index)
set(bwa.index)

len(.intersection())