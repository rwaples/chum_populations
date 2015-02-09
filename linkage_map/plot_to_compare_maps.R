library(ggplot2)



map_positions <- read.table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/test_union.txt", sep = "\t", header = TRUE)


offset = 73
map_positions$fam_08_LG <- map_positions$fam_08_LG-offset
map_positions$consensus_LG <- map_positions$consensus_LG-offset

ggplot(data = map_positions) + geom_point(aes(x = factor(fam_08_LG), y = fam_08_cM), color = 'blue', alpha = .5) +
  geom_point(aes(x = consensus_LG + .7, y = consensus_cM), color = 'red') +
  geom_segment(aes(x = fam_08_LG, y = fam_08_cM, xend =  consensus_LG + .7 , yend = consensus_cM), color = 'red', alpha = .2) +
  theme_bw()