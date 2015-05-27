library(ggplot2)

chum_ssa = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/synteny/chum_ssa.txt', 
                      sep ='\t', header = TRUE)

chum_ssa = chum_ssa[order(chum_ssa$LEP_LG),]
chum_ssa$chum_order = 1:nrow(chum_ssa)

#chum_ssa = chum_ssa[order(chum_ssa$Ssa_chr,chum_ssa$POS),]
chum_ssa = chum_ssa[order(chum_ssa$Ssa_chr),]
chum_ssa$ssa_order = 1:nrow(chum_ssa)

ggplot(chum_ssa) + geom_point(aes(x= chum_order, y = Ssa_chr, color = factor(LEP_LG%%5)), size = 7, alpha = .6) + 
  scale_x_continuous(name = 'chum LG', breaks = c(1,1+which(diff(chum_ssa$LEP_LG)!=0)), labels= 0:37) +
  scale_color_brewer(palette = 2, type = 'qual') + theme_bw()



