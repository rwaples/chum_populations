library(ggplot2)
library(RColorBrewer)

FDR =.05

#load base data
map_fst = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/pop_analysis/map_fst.txt', header = TRUE,sep = "\t")



map_fst$outlier = map_fst$qval < FDR
map_fst$LG = map_fst$paper1_LG


ggplot(map_fst, aes(x=cM, y=BS_Fst)) + 
  facet_grid(LG ~ .) +
  geom_point(aes(color=outlier), alpha = .7, size = 3) +
  scale_color_manual(values=c("darkgrey","red")) +
  stat_smooth(fill="blue", colour="darkblue", size=2, alpha = 0.1, method = 'loess', span=0.3, n = 10) +
  theme_bw() + theme(strip.text = element_blank()) + 
  theme(strip.text = element_text(size=20))
ggsave('/home/ipseg/Desktop/waples/chum_populations/paper/figures/supplemental/Fst_across_the_genome.png', height=48, width=12)

ggplot(map_fst, aes(x=cM, y=WC_Fst)) + 
  facet_grid(LG ~ .) +
  geom_point(aes(color=outlier), alpha = .7, size = 3) +
  scale_color_manual(values=c("darkgrey","red")) +
  stat_smooth(fill="blue", colour="darkblue", size=2, alpha = 0.1, method = 'loess', span=0.5, n = 20) +
  theme_bw() + theme(strip.text = element_blank()) + 
  theme(strip.text = element_text(size=20))
ggsave('/home/ipseg/Desktop/waples/chum_populations/paper/figures/supplemental/Fst_across_the_genome.png', height=48, width=12)


names

xx = map_fst[map_fst$LG ==10,]$cM
yy = map_fst[map_fst$LG ==10,]$BS_Fst
#yy[yy > 10000] = 5
scatter.smooth(xx, yy, span = .1, degree = 2, evaluation =100)

xx = map_fst[map_fst$LG ==8,]$cM
yy = map_fst[map_fst$LG ==8,]$BS_Fst
#yy[yy > 10000] = 5
scatter.smooth(xx, yy, span = .5, degree = 2, evaluation =100)
