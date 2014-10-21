args <- commandArgs(TRUE)

distance_matrix_file = args [1]
k = args [2]

distance_matrix = data.matrix(read.table(distance_matrix_file, header = FALSE))

mds = cmdscale(distance_matrix, k)

write.table()

