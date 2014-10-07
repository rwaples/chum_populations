library(ggplot2)


ind_stats_file = "C:/Users/IPSEG/Desktop/Waples/chum_populations/ind_seq_stats.tsv"
ind_stats <- read.table(file = ind_stats_file, sep = "\t", header = TRUE)

ggplot(data = ind_stats, aes(x= log10(nreads))) + 
  geom_histogram() 

ggplot(data = ind_stats, aes(x = pop, y= log10(nreads), fill = factor(pop) )) + 
  geom_violin(scale = 'width', trim = TRUE, adjust = 1) + coord_flip() 
 
ggplot(data = ind_stats, aes(x = pop, y= log10(nreads), fill = factor(pop) )) + 
  geom_violin(scale = 'width', trim = TRUE, adjust = 1) + coord_flip()

above = table(ind_stats[ind_stats$nreads > 750000,]$pop)
below = table(ind_stats[ind_stats$nreads < 750000,]$pop)

plot(above/(above+below))

