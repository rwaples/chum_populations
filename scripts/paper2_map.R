library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)

# load map
raw_map <-  read.table('/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/final/PS_chum_map_2015.txt', header = TRUE)
raw_map$paralog <- str_detect(raw_map[['resolved_locus']], '_')
raw_map$paralog <- str_detect(raw_map[['resolved_locus']], '_')

raw_map$paralog[raw_map$paralog== TRUE] <- "paralog"
raw_map$paralog[raw_map$paralog== FALSE] <- "non-paralog"

# load centromereic loci
cent_loci <- read.table('/home/ipseg/Desktop/waples/chum_populations/linkage_map/centromere_placement/cent_loci.txt', header = TRUE)
cent_loci_merged <- merge(cent_loci, raw_map, by = 'LEPname')

cent_loci_merged <- ddply(.data = cent_loci_merged, .variables = .(paper1_LG.y), mutate, top_cent = max(cM))
cent_loci_merged <- ddply(.data = cent_loci_merged, .variables = .(paper1_LG.y), mutate, bot_cent = min(cM))

cent_boundaries <- ddply(cent_loci_merged, 'paper1_LG.x', summarize, max = max(top_cent), min = min(bot_cent))


# Add column of max cM per LG
raw_map <- ddply(.data = raw_map, .variables = .(paper1_LG), mutate, end_cM = max(cM))
raw_map <- ddply(.data = raw_map, .variables = .(paper1_LG), mutate, low_cM = min(cM))
raw_map$cM <- raw_map$cM - raw_map$low_cM 

cM_size <- ddply(raw_map, 'paper1_LG', summarize, max = max(cM))
cM_size$start = 0

pl1 <- ggplot() + 
  # background lines
  geom_segment(data = cM_size, aes(x = paper1_LG, y = start , xend = paper1_LG, yend = max ), alpha = .2, size =1, color = 'black') +
  #centromeres
  geom_segment(data = cent_boundaries,  aes(x = paper1_LG.x, y = min, xend = paper1_LG.x, yend = max), 
               color = '#228B22',  alpha = 1, size = 8) +
  # loci
  geom_point(data = raw_map,     aes(x = paper1_LG, y = cM, color = paralog, shape = paralog), alpha = .4, size = 5 ) +
  # set colors
  scale_color_manual(values = c("#0000FF", "#fd3c06"), guide = guide_legend(title = NULL)) + 
  scale_shape_manual(values = c(16, 18), guide = guide_legend(title = NULL)) +   
  # legend
  guides(colour = guide_legend(ncol = 2, title = NULL, override.aes = list(alpha = 1))) +
  # x axis
  #scale_x_continuous(limits=c(1, 37)) +
  # style
  theme_minimal() +
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),         
        text = element_text(size=24),
        legend.justification=c(.5,0), legend.position=c(.5,.85)) +
  # labels
  xlab('Linkage Group') + 
  ylab('cM')  +
  ylim(0,150) + 
  scale_x_continuous(limits=c(1,37), breaks = seq(1,37,2) )

pl1

ggsave('/home/ipseg/Desktop/waples/chum_populations/paper/figures/linkage_map.png', width= 18, height = 8, dpi = 400 )

cent_boundaries$mycolor = 1:13
# centromeres
ggplot() + geom_segment(data = cent_boundaries,  aes(x = paper1_LG.x, y = min, xend = paper1_LG.x, yend = max, color = as.factor(mycolor)), 
             alpha = 1, size = 8) +   scale_fill_discrete(c('#228B22'), guide = guide_legend()) 
#geom_point(data = AB_loci_yes_PSV,     aes(x = LG_factor, y = cM), alpha = .8, size = 5, color = 'blue') +
#geom_point(data = AA_AB_loci_yes_PSV, aes(x = LG_factor , y= cM), shape = 18, alpha = .8, size = 5, col = 'orange') +
#geom_point(data = AA_BC_loci_yes_PSV, aes(x = LG_factor , y= cM), shape = 18, alpha = .8, size = 5, col = 'orange') +
#geom_point(data = AB_AC_loci_yes_PSV, aes(x = LG_factor , y= cM), shape = 18, alpha = .8, size = 5, col = 'orange') +
#geom_point(data = AB_CD_loci_yes_PSV, aes(x = LG_factor , y= cM), shape = 18, alpha = .8, size = 5, col = 'orange') +


sum(raw_map$paralog)
names(raw_map)
