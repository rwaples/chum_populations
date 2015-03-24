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
      
            

lepmap_files = glob.glob("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/*.map")
consensus_map = Map('consensus', 'lepmap', lepmap_files)                        
consensus_map.rename_markers("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_translation.txt", 2)

with open('/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap', 'w') as OUTFILE:
    OUTFILE.write('LG\tlocus\tcM\n')
    written_markers = []
    for LG, markers in consensus_map.markers_of_LG.items():
        for mark in markers:
            if mark not in written_markers:
                OUTFILE.write("{}\t{}\t{}\n".format(LG, mark,consensus_map.cM_of_marker[mark]))
                written_markers.append(mark)

#write_union('//home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap', consensus_map,consensus_map)



#fam08_map = Map('fam_08', 'mst', '/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_08_iter1_mstmap.map')
consensus_map = Map('consensus', 'generic', '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap')                        
fam08_map = Map('fam_08', 'generic', '/home/ipseg/Desktop/waples/chum_populations/linkage_map/Map2.txt')



# drop small LGs
for LG, markers in fam08_map.markers_of_LG.items():
    if len(markers) < 20:
        print(LG)
        fam08_map.drop_LG(LG)
        
for LG, markers in consensus_map.markers_of_LG.items():
    if len(markers) < 20:
        print(LG)
        consensus_map.drop_LG(LG)



# marker names used by lepmap are just positions, the markers must be renamed
# make rename file
# rename

# drop duplicated loci
consensus_map.drop_markers(*[xx for xx in consensus_map.LG_of_marker.keys() if "_" in xx])
fam08_map.drop_markers(*[xx for xx in fam08_map.LG_of_marker.keys() if "_x2" in xx]) # notice this just drops cases were both paralogs are mapped
# now rename these loci removing their "_x1"


with open('/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/08_x1_catIDs.txt', 'w') as OUTFILE:
    for loc in fam08_map.LG_of_marker.keys():
        OUTFILE.write("{}\t{}\n".format(loc.replace("_x1", ""), loc))

fam08_map.rename_markers("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/08_x1_catIDs.txt", 2)

compare_maps(fam08_map, consensus_map)

matches_1, maches_2, mismatches = align_maps(fam08_map,consensus_map, tolerance = 5)

consensus_map.drop_markers(*mismatches)
fam08_map.drop_markers(*mismatches)

#consensus_map.drop_markers('3619')
#consensus_map.drop_markers('48610_NA')

align_maps(fam08_map,consensus_map, tolerance = 5)
rhos = compare_orders(fam08_map, consensus_map)
write_union('/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/partial_union.txt', fam08_map,consensus_map)

################################################

consensus_map.LG_of_marker.keys()  ['48610_NA']
consensus_map.LG_of_marker['48610_NA']
consensus_map.cM_of_marker['48610_NA']
consensus_map.cM_of_marker.keys()

align_maps(consensus_map,fam08_map)

rhos = compare_orders(consensus_map, fam08_map)


for loc in consensus_map.LG_of_marker.keys():
    Map.drop_markers([xx for xx in ])




consensus_map.LG_of_marker.keys()
fam08_map.LG_of_marker.keys()





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