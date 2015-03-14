library(ggplot2)


inital <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes'
lod7 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod7_singles.chromosomes'
lod6 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod6_singles.chromosomes'
lod5 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod5_singles.chromosomes'
lod4 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod4_singles.chromosomes'
lod3.5 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod3-5_singles.chromosomes'

collapsed <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.inital.chromosomes'
lod7 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod7_singles.chromosomes'
lod6 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod6_singles.chromosomes'
lod5 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes'
lod4 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod4_singles.chromosomes'
#lod3.5 <- '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod3-5_singles.chromosomes'


plot_LG_hist <- function(path_to_LEPmap_chr_file){
  chr_assigments = read.table(path_to_LEPmap_chr_file, header = FALSE)
  assigned_to_LG <- subset(chr_assigments, V1 > 0)
  ggplot(data = assigned_to_LG, aes(x = V1)) + geom_histogram(binwidth = 1, color = 'darkgreen', fill = 'grey') + 
    geom_text(label = paste("Unassigned: ", as.character(nrow(chr_assigments)-nrow(assigned_to_LG))), x = 30, y = max(table(assigned_to_LG$V1))-40) + theme_minimal()
} 

plot_LG_hist(inital)
plot_LG_hist(collapsed)
plot_LG_hist(lod7)
plot_LG_hist(lod6)
plot_LG_hist(lod5)
plot_LG_hist(lod4)
plot_LG_hist(lod3.5)

# LG congruence
library(ggplot2)
library(igraph)
library(reshape2)
LG_congru <- read.table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/LG_congruence.tsv", header = FALSE, sep = "\t")
names(LG_congru) <- c('locus', 'family', 'copy', 'LG')
LG_congru = LG_congru[LG_congru$LG != 0,]

rename_table <- read.table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/rename_table.tsv", header = FALSE, sep = "\t")
names(rename_table) <- c('fam_name', 'name', 'LG')
rename_table = rename_table[rename_table$LG != 0,]

adj <- crossprod(table(rename_table[c(2,3)]))
#diag(adj) <- 0
sum(adj)
adj_matrix_melt <- melt(adj)
names(adj_matrix_melt) <- c('LG_1', 'LG_2', 'overlap')
adj_matrix_melt <- adj_matrix_melt[adj_matrix_melt$overlap > 1,]
ggplot(data = adj_matrix_melt, aes(x = factor(LG_1), y = factor(LG_2))) + geom_tile(aes(fill = overlap)) + 
  scale_fill_continuous(low="yellow", high="red")


occur_multiple = names(table(LG_congru$locus) [table(LG_congru$locus) > 1])
LG_multiple_congru = LG_congru[LG_congru$locus %in% occur_multiple, ]


adj <- crossprod(table(LG_congru[c(1,4)]))
diag(adj) <- 0
g <- graph.adjacency(adj, weighted=TRUE, mode ='undirected', diag = FALSE)
g <- graph.adjacency(adj, weighted=TRUE, mode ='min', diag = FALSE)
g <- simplify(g)
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
plot(g)

adj_matrix_melt <- melt(adj)
names(adj_matrix_melt) <- c('LG_1', 'LG_2', 'overlap')
adj_matrix_melt <- adj_matrix_melt[adj_matrix_melt$overlap > 1,]
ggplot(data = adj_matrix_melt, aes(x = factor(LG_1), y = factor(LG_2))) + geom_tile(aes(fill = overlap)) + 
  scale_fill_continuous(low="yellow", high="red")

heatmap((adj[upper.tri(adj)]))
adj[upper.tri(adj)]
data.frame(adj_matrix_melt[upper.tri(adj_matrix_melt)])


adj

bip <- graph.data.frame(LG_multiple_congru[c("LG", "locus")])
plot(bip)
V(bip)$type <- V(bip)$name %in% as.character(LG_multiple_congru[,1])
v <- get.adjacency(bipartite.projection(bip)[[2]], attr="weight", sparse=FALSE)


library(reshape2)
melt(table(LG_multiple_congru$locus, LG_multiple_congru$LG))

aa = melt(LG_multiple_congru[c('locus', 'LG')], id = 'locus')

ggplot(data = aa, aes(x = locus)) + geom_point()

heatmap(cov(table(LG_multiple_congru$locus, LG_multiple_congru$LG)))

table(LG_congru$locus, LG_congru$LG)

duplicated(LG_congru$locus)[1:10]
LG_congru$locus[1:10]



heatmap(cor(table(LG_multiple_congru$LG, LG_multiple_congru$locus)))


(LG_congru)
