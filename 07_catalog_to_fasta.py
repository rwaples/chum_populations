# convert a stacks catalog to a reference fasta 


""" 
TODO: filter catalog based on:
    repetitive sequences
    unique alignments
    rxstacks
    mappability
    depth



"""




#my_tags_file = '/media/Shared/Data/chum/populations/stacks/batch_1.catalog.tags.tsv'
my_ref_file = '/media/Shared/Data/chum/populations/aln/chum_ref_batch_03.fa'
my_tags_file = "/media/Shared/Data/chum/populations/stacks/batch_3.catalog.tags.tsv"
#my_ref_file = '/media/Shared/Data/chum/populations/aln/chum_ref_batch_02.fa'

def tags_to_ref(tags_file, ref_file):
    with open(tags_file) as INFILE:
        with open(ref_file, 'w') as OUTFILE: 
            for line in INFILE:
                catalog_split = line.strip().split()
                catID = catalog_split[2]
                seq = catalog_split[8]
                OUTFILE.write(">catID|{}\n{}\n".format(catID, seq))
                
tags_to_ref(my_tags_file, my_ref_file)

