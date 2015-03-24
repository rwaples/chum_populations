library(ggplot2)
library(stringr)
library(plyr)


map_positions <- read.table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/partial_union.txt", sep = "\t", header = TRUE)


# offset = 73
# map_positions$fam_08_LG <- map_positions$fam_08_LG-offset
# map_positions$consensus_LG <- map_positions$consensus_LG-offset

ggplot(data = map_positions) + geom_point(aes(x = factor(fam_08_LG), y = fam_08_cM), color = 'blue', alpha = .5) +
  geom_point(aes(x = consensus_LG +.2, y = consensus_cM), color = 'red') +
  geom_segment(aes(x = fam_08_LG, y = fam_08_cM, xend =  consensus_LG + .2 , yend = consensus_cM), color = 'red', alpha = .2) +
  theme_bw()


collapsed_map_positions <- read.table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap", sep = "\t", header = TRUE)
names(collapsed_map_positions)

split_character<- str_split_fixed(collapsed_map_positions$Marker, "_", 2)[,2]
collapsed_map_positions$duplicated  = (split_character == "A") | (split_character == "B")
collapsed_map_positions <- collapsed_map_positions[collapsed_map_positions$consensus_LG <= 37,]
collapsed_map_positions$catID <- str_split_fixed(collapsed_map_positions$Marker, "_", 2)[,1]
lenLG <- ddply(collapsed_map_positions, 'consensus_LG', summarize, lenLG = max(consensus_cM))

paralogs = collapsed_map_positions[collapsed_map_positions$duplicated == TRUE,]
paralogs$catID <-  str_split_fixed(paralogs$Marker, "_", 2)[,1]
homeologs <- merge(paralogs, paralogs, by = 'catID', all = FALSE, stringsAsFactors = FALSE)
homeologs <- homeologs[homeologs$Marker.x != homeologs$Marker.y,]
names(homeologs)
#mi <- homeologs[seq(1, nrow(homeologs), 2), ]


ggplot(data = collapsed_map_positions) + 
  geom_segment(data = lenLG, aes(x = factor(consensus_LG), xend = factor(consensus_LG), y = 0, yend = lenLG), alpha = .2, size =1, color = 'black') +
  geom_point(aes(x = factor(consensus_LG), y = consensus_cM, color = duplicated, shape=duplicated, order =duplicated, size = duplicated), alpha = .5) + 
  scale_colour_manual(values = c("blue","orange")) + scale_shape_manual(values=c(16,18)) + scale_size_manual(values = c(5,6)) +
  geom_segment(data = homeologs, aes(x = consensus_LG.x, xend = consensus_LG.y, y = consensus_cM.x, yend=consensus_cM.y) ) +
  theme_minimal() +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        text = element_text(size=24))

write.table(unique(collapsed_map_positions$catID), file = "/home/ipseg/Desktop/waples/chum_populations/data/batch_2/on_map.txt", quote = FALSE, row.names = FALSE)


sum(collapsed_map_positions$duplicated)
sum(homeologs$duplicated.y)/4
