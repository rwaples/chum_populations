library(ggplot2)


ind_stats_file = file.path('./results/ind_seq_stats.tsv')
ind_stats <- read.table(file = ind_stats_file, sep = "\t", header = TRUE)

ind_stats = data.frame(ind_stats[ind_stats$pop != 'CMUW10X', ])
ind_stats = data.frame(ind_stats[ind_stats$pop != 'MUWhap10', ])
ind_stats$pop = factor(ind_stats$pop) # reset levels
nrows = nrow(ind_stats)

# histogram of depth
ggplot(data = ind_stats, aes(x= log10(nreads))) + 
  geom_histogram(binwidth = .2) + geom_vline(xintercept = 6, color = 'red') + theme_classic()

# indidivuals ordered by depth
ggplot(data = ind_stats[order(ind_stats$nreads),], aes(y= log10(nreads), x = seq(1, nrows))) +
  geom_point(size = 5, alpha = .8) + geom_hline(yintercept = 6, color = 'black') + 
  xlab('individual') +
  theme_classic() + theme(axis.text = element_text(size = 14, face = 'bold'), axis.title = element_text(size = 24, face = 'bold')) 

# depth by population
ggplot(data = ind_stats, aes(x = pop, y= log10(nreads), fill = factor(pop) )) + 
  geom_violin(scale = 'width', trim = TRUE, adjust = 1) + geom_hline(yintercept = 6, color = 'black') + coord_flip() + theme_classic() + theme(axis.text = element_text(size = 14, face = 'bold'), axis.title = element_text(size = 24, face = 'bold')) 
 
above = table(ind_stats[ind_stats$nreads > 1000000,]$pop)
below = table(ind_stats[ind_stats$nreads < 1000000,]$pop)
above_below = data.frame(cbind(above, below))
above_below$pop = rownames(above_below)

# sample size at depth thresholds
ggplot(data = above_below, aes(x = below, y = above, fill = factor(pop), color = factor(pop))) + 
  geom_point(size = 8, alpha = .5) + geom_text(aes(label = pop), hjust=0, vjust=0, size = 10, position = 'jitter') + 
  xlab("# individuals below 1M depth") + ylab("# individuals above 1M depth") + 
  theme_classic() + xlim(-1, 13) + theme(axis.text = element_text(size = 14, face = 'bold'), axis.title = element_text(size = 24, face = 'bold')) 

plot(log10(sort(ind_stats$nreads)))




plot(above/(above+below))


plot(as.list(above), as.list(below), xlab = 'individuals above 1mil depth', type = 'n')
text(as.list(above), as.list(below), names(above))


levels(ind_stats$pop)
