library(ggplot2)
library(RColorBrewer)
args <- commandArgs(TRUE)

mds_file = args[1]
mds = read.table(mds_file, sep = "\t", header = TRUE)
mds$TIMING <- factor(mds$TIMING, levels = c('Summer', "Fall", "Winter", "UNK"))
#mds = read.table("./data/batch_10/batch_10.dist_mds", sep = "\t", header = TRUE)

my_cols <- c("#000000", brewer.pal(9, "Set1"))

plot1v2 <- ggplot(data = mds, aes(x = X1, y = X2)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(POPNAME), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  theme_bw() + coord_equal() + scale_color_manual(values = my_cols) +
  xlab("Axis 1 - 2.76%") + ylab("Axis 2 - 2.47%") +
  ggtitle("Puget Sound Chum Salmon Populations")

ggsave(filename = "./plots/1v2.pdf", plot1v2, width = 12, height = 12)

plot1v2region <- ggplot(data = mds, aes(x = X1, y = X2)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(REGION), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  #  geom_segment(data= batch_10_weights, mapping = aes(x = 0, y= 0, xend = Dim1/3, yend = Dim2/3)) +
  theme_bw() + coord_equal() + scale_color_brewer(type = 'qual', palette="Set1") +
  xlab("Axis 1 - 2.76%") + ylab("Axis 2 - 2.47%") +
  ggtitle("Puget Sound Chum Salmon Populations")
ggsave(filename = "./plots/1v2region.pdf", plot1v2region, width = 12, height = 12)

plot1v3 <-ggplot(data = mds, aes(x = X1, y = X3)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(POPNAME), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  theme_bw() + coord_equal() + scale_color_manual(values = my_cols) +
  xlab("Axis 1 - 2.76%") + ylab("Axis 3 - 1.69%") +
  ggtitle("Puget Sound Chum Salmon Populations")
ggsave(filename = "./plots/1v3.pdf", plot1v3, width = 12, height = 12)

plot1v3region <-ggplot(data = mds, aes(x = X1, y = X3)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(REGION), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  #  geom_segment(data= batch_10_weights, mapping = aes(x = 0, y= 0, xend = Dim1/3, yend = Dim3/3)) +
  theme_bw() + coord_equal() +  scale_color_brewer(type = 'qual', palette="Set1") +
  xlab("Axis 1 - 2.76%") + ylab("Axis 3 - 1.69%") +
  ggtitle("Puget Sound Chum Salmon Populations")
ggsave(filename = "./plots/1v3region.pdf", plot1v3region, width = 12, height = 12)

plot2v3 <-ggplot(data = mds, aes(x = X2, y = X3)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(POPNAME), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  theme_bw() + coord_equal() + scale_color_manual(values = my_cols) +
  xlab("Axis 2 - 2.47%") + ylab("Axis 3 - 1.69%") +
  ggtitle("Puget Sound Chum Salmon Populations")
ggsave(filename = "./plots/2v3.pdf", plot2v3, width = 12, height = 12)

plot2v3region <-ggplot(data = mds, aes(x = X2, y = X3)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 8, alpha = .8, mapping = aes(color = factor(REGION), shape = factor(TIMING))) +
  geom_point(size = 3, alpha = .8, mapping = aes(shape = factor(TIMING))) +
  theme_bw() + coord_equal() +  scale_color_brewer(type = 'qual', palette="Set1") +
  xlab("Axis 2 - 2.47%") + ylab("Axis 3 - 1.69%") +
  ggtitle("Puget Sound Chum Salmon Populations")
ggsave(filename = "./plots/2v3region.pdf", plot2v3region, width = 12, height = 12)
