import glob
from parseMap import *

with open("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap")as INFILE:
    #just need first line
    first_line = next(INFILE)
    marker_names = first_line.strip().split()[6:]
    
with open("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_translation.txt", 'w') as OUTFILE:
    count = 0
    for marker in marker_names:
        count += 1
        OUTFILE.write( "{}\t{}\n".format(marker, count))
            
            


fam08_map = Map('fam_08', 'mst', '/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_08_iter1_mstmap.map')

lepmap_files = glob.glob("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/*.map")
consensus_map = Map('consensus', 'lepmap', lepmap_files)


# marker names used by lepmap are just positions, the markers must be renamed
# make rename file
# rename
consensus_map.rename_markers("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_translation.txt", 2)
write_union('//home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap', consensus_map,consensus_map)


consensus_map.LG_of_marker.keys()
fam08_map.LG_of_marker.keys()
compare_maps(fam08_map, consensus_map)




fam08_map.drop_LG('1')
consensus_map.drop_LG('38')
consensus_map.drop_LG('39')
consensus_map.drop_LG('40')
consensus_map.drop_LG('41')

align_maps(fam08_map,consensus_map)
align_maps(consensus_map,fam08_map)

compare_orders(fam08_map, consensus_map)
rhos = compare_orders(consensus_map, fam08_map)

for LG in fam08_map.markers_of_LG.keys():
    if LG.startswith("un_"):
        fam08_map.drop_LG(LG)
        
for LG in consensus_map.markers_of_LG.keys():
    if LG.startswith("un_"):
        consensus_map.drop_LG(LG)


write_union('/home/ipseg/Desktop/waples/chum_populations/linkage_map/test_union.txt', fam08_map,consensus_map)
write_union('/home/ipseg/Desktop/waples/chum_populations/linkage_map/centromere_placement/fam_08_map.txt', fam08_map,fam08_map)


rhos = compare_orders(consensus_map, fam08_map)
len([x for x,y in rhos.items() if y > .90])