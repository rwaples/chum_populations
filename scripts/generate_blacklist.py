import pandas

markers_filename = '/media/Shared/Data/chum/populations/stacks_output/populations/batch_10.markers.tsv'

markers = pandas.DataFrame.from_csv(markers_filename, sep = "\t")



genotyped_loci = markers.loc[markers["Total Genotypes"] > 0,]["Catalog Locus ID"]
genotyped_loci.to_csv('/media/Shared/Data/chum/populations/stacks_output/populations/whitelist.txt')



markers["Total Genotypes"][0,]