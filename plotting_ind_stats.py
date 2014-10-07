library(ggplot2)


ind_stats_file = "C:/Users/IPSEG/Desktop/Waples/chum_populations/ind_seq_stats.tsv"
ind_stats <- read.table(file = ind_stats_file, sep = "\t", header = TRUE)


ggplot(data = ind_stats, aes(x = pop, y= log10(nreads))) + 
  geom_boxplot() 
  
  fill = pop, color = pop   
  
ggplot(data = ind_stats, aes(x = pop, y= lenreads, fill = pop, color = pop )) + 
  geom_boxplot()
  
  facet_grid(popsize_factor ~ starting_loci_factor, labeller = label_both, scales = "free_y") + 
  labs(title = "Effecive Size estimated from all loci", y= 'Ne', x= "Number of chromosomes")