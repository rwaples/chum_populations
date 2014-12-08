#evecs  <-  read.table("C:/Users/Kele/Desktop/chum_populations/data/test_EIGENSTRAT/FISH_560.genotypes.evec", header = FALSE, skip = 1)
library(stringr)
library(ggplot2)
library(candisc)

eigenvec <- read.table("C:/Users/Kele/Desktop/chum_populations/data/FISH_560/FISH_560.final.eigenvec", sep = "\t", header = TRUE)


eigenvec$POP <- str_split_fixed(eigenvec[,2], "_", 2)[,1]
popdata = read.table("C:/Users/Kele/Desktop/chum_populations/data/POPINFO.txt", sep = "\t", header = TRUE)
names(popdata)[1] <- 'POP'

eigenvec <- merge(eigenvec, popdata, all = TRUE )
eigenvec$TIMING <- factor(eigenvec$TIMING, levels = c('Summer', "Fall", "Winter"))

# try kicking out winters
#eigenvec <- eigenvec[eigenvec$TIMING != 'Winter',]

my.mod.REGION <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ (REGION + TIMING) , data=eigenvec)
my.can.REGION <- candisc(my.mod.REGION, data=eigenvec, term = "REGION")
names(my.can.REGION)

my.mod.TIMING <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ (REGION + TIMING), data=eigenvec)
my.can.TIMING <- candisc(my.mod.TIMING, data=eigenvec, "TIMING")
names(my.mod.TIMING)

my.mod.POP <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ POP, data=eigenvec)
my.can.POP <- candisc(my.mod.POP, data=eigenvec)
names(my.can.POP)


theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}

CDA_REGION <- ggplot() + geom_point(data =  my.can.REGION$scores, aes(x = Can1, y = Can2, color = REGION, shape = TIMING), size = 3 , alpha = .8) + 
  geom_point(data =  my.can.REGION$scores, aes(x = Can1, y = Can2, shape = TIMING), color = 'black', size = 1.5 , alpha = .5) +
  geom_segment(data = data.frame(my.can.REGION$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 1.5, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.REGION$coeffs.std, row.names(my.can.REGION$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.REGION.coeffs.std.), size = 4) + 
  xlab("Can1 - 61%") + ylab("Can2 - 38%") + coord_equal() + theme_bw() 

ggsave(filename = "C:/Users/Kele/Desktop/chum_populations/plots/CDA_REGION.pdf", CDA_REGION, width = 6, height = 6)

CDA_TIMING <- ggplot() + geom_point(data =  my.can.TIMING$scores, aes(x = Can1, y = Can2, color = REGION, shape = TIMING), size = 3 , alpha = .8) + 
  geom_point(data =  my.can.TIMING$scores, aes(x = Can1, y = Can2, shape = TIMING), color = 'black', size = 1.5 , alpha = .5) +
  geom_segment(data = data.frame(my.can.TIMING$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 1.5, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.TIMING$coeffs.std, row.names(my.can.TIMING$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.TIMING.coeffs.std.), size = 4) + 
  xlab("Can1 - 92%") + ylab("Can2 - 7.4%") +
  theme_bw()

ggsave(filename = "C:/Users/Kele/Desktop/chum_populations/plots/CDA_TIMING.pdf", CDA_TIMING, width = 6, height = 6)

ggplot() + geom_point(data =  my.can.POP$scores, aes(x = Can1, y = Can2, color = POP), size = 8 , alpha = .8) + 
  geom_segment(data = data.frame(my.can.POP$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 2, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.POP$coeffs.std, row.names(my.can.POP$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.POP.coeffs.std.), size = 10) + 
#  xlab("Can1 - 61%") + ylab("Can2 - 38%") +
  theme_bw()


# test
my.mod.REGION <- lm(cbind( PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8, PC9) ~ (POP) , data=eigenvec)
my.can.REGION <- candisc(my.mod.REGION, data=eigenvec, term = "POP")
plot(my.can.TIMING, suffix = TRUE)

ggplot() + geom_point(data =  my.can.REGION$scores, aes(x = Can1, y = Can2, color = REGION, shape = TIMING), size = 8 , alpha = .8) + 
  geom_segment(data = data.frame(my.can.REGION$coeffs.std), x = 0, y = 0, aes(xend = Can1*5 , yend = Can2*5), size = 2, alpha = .4,color = "blue") +
  geom_text(data = data.frame(my.can.REGION$coeffs.std, row.names(my.can.REGION$coeffs.std)), aes(x = Can1*5 , y = Can2*5, label = row.names.my.can.REGION.coeffs.std.), size = 10) + 
  xlab("Can1 - 61%") + ylab("Can2 - 38%") +
  theme_bw()

library(vegan)
cc = eigenvec[-142,]
aa = with(cc, cbind(PC1,PC2,PC3,PC4,PC5,PC6))
bb = with(cc, cbind(POP,REGION, TIMING))
cca(aa, bb)
  
    )