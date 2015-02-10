# place centromeres
import glob
from collections import OrderedDict

from parseMap import *





fam08_map = Map('fam_08', 'mst', '/home/ipseg/Desktop/waples/chum_populations/linkage_map/centromere_placement/chum_08_mstmap.map')

def parse_genotypes(genotypes_file):
    with open(genotypes_file) as INFILE:
        ind_names = next(INFILE).strip().split("\t")[1:]
        genotypes_at_locus = dict()
        for line in INFILE:
            locus = line.strip().split("\t")[0]
            genotypes = line.strip().split("\t")[1:]
            genotype_of_ind = OrderedDict(zip(ind_names, genotypes))
            genotypes_at_locus[locus] = genotype_of_ind
        
    return(genotypes_at_locus)
    
# genotypes format - rows: loci columns: individuals
fam08_genotypes_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/centromere_placement/chum_08_genotypes.txt" 
genos = parse_genotypes(fam08_genotypes_file)


LG_order = sorted(fam08_map.markers_of_LG.keys(), key = int)
marker_order = list()
for LG in LG_order:
    if len(fam08_map.markers_of_LG[LG]) > 1:
        markers_on_LG = sorted(fam08_map.markers_of_LG[LG], key=lambda k: fam08_map.cM_of_marker[k])
        marker_order.extend(markers_on_LG)
        #marker_order.extend(sorted(fam08_map.markers_of_LG[LG], key=lambda k: fam08_map.cM_of_marker[k]))       


def write_centro_place_file(input_Map, input_genos, filename):
    # establish marker order
    LG_order = sorted(input_Map.markers_of_LG.keys(), key = int)
    marker_order = list()
    for LG in LG_order:
        if len(input_Map.markers_of_LG[LG]) <= 2:
            print("excluding LG: {}".format(LG))
        else:
            markers_on_LG = input_Map.markers_of_LG[LG]
            cM_order = sorted(markers_on_LG, key=lambda k: input_Map.cM_of_marker[k])
            marker_order.extend(cM_order)       
    with open(filename, 'w') as OUTFILE:
        # write locus names
        OUTFILE.write("locus_names," +",".join(marker_order) + "\n")
        OUTFILE.write("LG," +",".join([str(input_Map.LG_of_marker[xx]) for xx in marker_order]) + "\n")
        OUTFILE.write("cM," +",".join([str(input_Map.cM_of_marker[xx]) for xx in marker_order]) + "\n")
        OUTFILE.write("cM_centro," +",".join(["NA" for xx in marker_order]) + "\n")
        OUTFILE.write("Chr_Arm," +",".join(["NA" for xx in marker_order]) + "\n")
        OUTFILE.write("Y," +",".join(["NA" for xx in marker_order]) + "\n")
        for ind in genos.values()[0].keys():
            OUTFILE.write("{},".format(ind) + ",".join([genos[xx][ind] for xx in marker_order]) + "\n")




write_centro_place_file(fam08_map, genos, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/centromere_placement/chum_08_input.csv" )
