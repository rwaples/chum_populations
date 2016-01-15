source("https://bioconductor.org/biocLite.R")
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

biocLite("ggtree")

library("ggtree")

nwk_file = "/home/ipseg/Desktop/waples/chum_populations/results/batch_4/pop_analysis/poptree2/FST_named.nwk"
tree <- read.tree(nwk_file)

tr <- ggtree(tree, layout="unrooted") + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
tr
tr %>% rotate(11)

tr <- ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
tr
flip(tr, 17, 12 ) 


tr <- ggtree(tree, layout="unrooted") + ggtitle("unrooted layout") + geom_label(aes(label=node)) 
tr
flip(tr, 17, 11 ) 







# Find the internal node number of the input taxas MRCA
MRCA(tree, tip = c('Hoodsport_Hatchery', 'Lilliwaup_Creek'))


ggtree(tree, layout="unrooted") + ggtitle("unrooted layout")    +
#geom_tiplab(size=6, color="black")      +
  geom_label(aes(label=node))   %>% flip(12)  +
#  geom_text(aes(label=label))  
#  geom_label(aes(label=label))     +
  add_legend()  + 
  theme_tree()

ggtree(tree, branch.length="none")



%>% rotate(21)



geom_tiplab(aes(angle=angle), size=3, color="black")
+ geom_treescale(fontsize=8, linesize=2, offset=0)
