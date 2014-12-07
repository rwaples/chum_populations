library(stringr)
args <- commandArgs(TRUE)

POPINFO_file = args[1]
ind_ID_file = args[2]
distance_matrix_file = args[3]
k = args[4]
file_out = args[5]
file_out = paste(args[5], ".eigv", sep = "")

#popinfo = read.table("POPINFO.txt", sep = "\t", header = TRUE)
#ind_IDs = read.table("./data/batch_10/batch_10.mdist.id", sep = "\t", header = FALSE)

popinfo = read.table(POPINFO_file, sep = "\t", header = TRUE)
ind_IDs = read.table(ind_ID_file, sep = "\t", header = FALSE)

names(ind_IDs) <- c('POPNUM', 'IND')
ind_IDs$SILLI <- str_split_fixed(ind_IDs$IND, '_', 2)[,1]
ind_IDs_with_POPINFO <- merge(ind_IDs, popinfo)


distance_matrix = data.matrix(read.table(distance_matrix_file, header = FALSE))
mds = data.frame(cmdscale(distance_matrix, k))
mds$POPNAME <- ind_IDs_with_POPINFO$POPNAME
mds$INS <- ind_IDs_with_POPINFO$IND
mds$TIMING <- ind_IDs_with_POPINFO$TIMING
mds$REGION <- ind_IDs_with_POPINFO$REGION
  
write.table(mds, file = file_out, sep = "\t", quote = FALSE, row.names = FALSE)

eigenvals <- cmdscale(distance_matrix, k = k, eig = TRUE)$eig

write.table(eigenvals, file = file_out, sep = "\t", quote = FALSE, row.names = FALSE)

