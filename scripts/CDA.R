#evecs  <-  read.table("C:/Users/Kele/Desktop/chum_populations/data/test_EIGENSTRAT/FISH_560.genotypes.evec", header = FALSE, skip = 1)

eigenvec <- read.table("C:/Users/Kele/Desktop/chum_populations/data/FISH_560/FISH_560.final.eigenvec", sep = "\t", header = TRUE)

library(stringr)
eigenvec$pop <- str_split_fixed(eigenvec[,2], "_", 2)[,1]
popdata = read.table("C:/Users/Kele/Desktop/chum_populations/data/POPINFO.txt", sep = "\t", header = TRUE)
names(popdata)[1] <- 'pop'

eigenvec <- merge(eigenvec, popdata, all = TRUE )

my.mod.REGION <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ REGION, data=eigenvec)
my.can.REGION <- candisc(my.mod.REGION, data=eigenvec)
names(my.can.REGION)

my.mod.TIMING <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ TIMING, data=eigenvec)
my.can.TIMING <- candisc(my.mod.TIMING, data=eigenvec)
names(my.mod.TIMING)

library(ggplot2)
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}

ggplot() + geom_point(data =  my.can.REGION$scores, aes(x = Can1, y = Can2, color = REGION), size = 8 , alpha = .8) + 
  geom_segment(data = data.frame(my.can.REGION$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 2, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.REGION$coeffs.std, row.names(my.can$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.REGION.coeffs.std.), size = 10) + 
  xlab("Can1 - 61%") + ylab("Can2 - 38%") +
  theme_bw()


ggplot() + geom_point(data =  my.can.TIMING$scores, aes(x = Can1, y = Can2, color = TIMING), size = 8 , alpha = .8) + 
  geom_segment(data = data.frame(my.can.TIMING$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 2, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.TIMING$coeffs.std, row.names(my.can.TIMING$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.TIMING.coeffs.std.), size = 10) + 
  xlab("Can1 - 61%") + ylab("Can2 - 38%") +
  theme_bw()



