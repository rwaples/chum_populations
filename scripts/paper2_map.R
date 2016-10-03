library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)

# load map
raw_map <-  read.table('C:/Users/Ryan/Documents/GitHub/chum_populations/linkage_map/LEPmap/with_paralogs/final/PS_chum_map_2015.txt', header = TRUE)
raw_map$paralog <- str_detect(raw_map[['resolved_locus']], '_')
raw_map$paralog <- str_detect(raw_map[['resolved_locus']], '_')

raw_map$paralog[raw_map$paralog== TRUE] <- "paralog"
raw_map$paralog[raw_map$paralog== FALSE] <- "non-paralog"

# load centromereic loci
cent_loci <- read.table('C:/Users/Ryan/Documents/GitHub/chum_populations/linkage_map/centromere_placement/cent_loci.txt', header = TRUE)
cent_loci_merged <- merge(cent_loci, raw_map, by = 'LEPname')

cent_loci_merged <- ddply(.data = cent_loci_merged, .variables = .(paper1_LG.y), mutate, top_cent = max(cM))
cent_loci_merged <- ddply(.data = cent_loci_merged, .variables = .(paper1_LG.y), mutate, bot_cent = min(cM))

cent_boundaries <- ddply(cent_loci_merged, 'paper1_LG.x', summarize, max = max(top_cent), min = min(bot_cent))
# make sure they are at least 1 cm wide
cent_boundaries[cent_boundaries$min == cent_boundaries$max,]$max = cent_boundaries[cent_boundaries$min == cent_boundaries$max,]$min + 1


# Add column of max cM per LG
raw_map <- ddply(.data = raw_map, .variables = .(paper1_LG), mutate, end_cM = max(cM))
raw_map <- ddply(.data = raw_map, .variables = .(paper1_LG), mutate, low_cM = min(cM))
raw_map$cM <- raw_map$cM - raw_map$low_cM 

cM_size <- ddply(raw_map, 'paper1_LG', summarize, max = max(cM))
cM_size$start = 0


# Supplemental figure
sup1 <- ggplot() + 
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

sup1

#ggsave('C:/Users/Ryan/Documents/GitHub/chum_populations/paper/figures/linkage_map.png', width= 16, height = 8, dpi = 400 )



### Main figure 2 ###
to_plot = c(16, 29)

to_plot_cM_size = cM_size[cM_size[['paper1_LG']] %in% to_plot,]
to_plot_cent_boundaries = cent_boundaries[cent_boundaries[['paper1_LG.x']] %in% to_plot,]
to_plot_raw_map = raw_map[raw_map[['paper1_LG']] %in% to_plot,]
to_plot_raw_map['cM_rev'] = to_plot_raw_map['end_cM'] - to_plot_raw_map['cM']


testmerge = merge(to_plot_raw_map, to_plot_raw_map, by = c('contig'))
testmerge1 = testmerge[testmerge[['paper1_LG.x']] != testmerge[['paper1_LG.y']],]
testmerge2 = testmerge1[!duplicated(testmerge1['resolved_locus.x']),]

main2a <- ggplot() + 
  # background lines
  geom_segment(data = to_plot_cM_size, aes(y = paper1_LG, x = start , yend = paper1_LG, xend = max ), alpha = .2, size =1, color = 'black') +
  #centromeres
  geom_segment(data = to_plot_cent_boundaries,  aes(y = paper1_LG.x, x = min, yend = paper1_LG.x, xend = max), 
               color = '#228B22',  alpha = 1, size = 8) +
  # paralog connections
  geom_segment(data = testmerge2, aes(y = paper1_LG.x, x = cM.x , yend = paper1_LG.y, xend = cM.y ), alpha = .1, size =1, color = 'black') +
  
# loci
  geom_point(data = to_plot_raw_map, aes(y = paper1_LG, x = cM, color = paralog, shape = paralog), alpha = .4, size = 5 ) 

# set colors
main2b = main2a +scale_color_manual(values = c("#0000FF", "#fd3c06"), 
                     guide = guide_legend(title = NULL)
                     ) + 
  scale_shape_manual(values = c(16, 18), 
                     guide = guide_legend(title = NULL)
                     ) +   
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
        legend.justification=c(.5,0), 
        legend.position=c(.5,.85),
        plot.margin = margin(20, 10, 20, 30),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text.y=element_text(margin=margin(0,-20,0,20)),
        axis.text.x=element_text(margin=margin(-10,0,0,20))        
        ) +
  # labels
  ylab('Linkage\nGroup') + 
  xlab('cM')  +
  scale_y_continuous(limits=c(12,40), breaks = c(16,29)) + 
  scale_x_continuous(limits=c(0,130), breaks = c(0, 25, 50, 75, 100, 125))

main2b
ggsave('C:/Users/Ryan/Documents/GitHub/chum_populations/paper/figures/figure2.png', width= 6, height = 4, dpi = 400 )

##### END ####





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
