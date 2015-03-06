from multi_family_mapping_functions import *
import pandas as pd

# Settings
#os.chdir("Y:/WORK/WAPLES/Stacks_mapping/Python")

linkage_map_file_1 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_08_mstmap.txt"
linkage_map_file_2 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_01_mstmap.txt"
linkage_map_file_3 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_09_mstmap.txt"

stats_file_1 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_08.stats"
stats_file_2 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_01.stats"
stats_file_3 = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_09.stats"
blacklist_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_paralogs.txt"

   
# create a file listing duplicated loci (blacklist)
# these are loci, that are listed as confounded in any of the three familes 
for stats_file in [stats_file_1, stats_file_2, stats_file_3]:
    paralogs = set()
    with open(stats_file) as INFILE:
        for line in INFILE:
            if line.strip().split()[2] not in ['AA_xx', "AB"]:
                paralogs.add(line.strip().split()[0])
    with open(blacklist_file, 'w') as OUTFILE:
        for xx in paralogs:
            OUTFILE.write(xx)
            OUTFILE.write("\n")
                                                                                                      
# Import genotype data from MST map input files:
individuals_08, genotypes_at_locus_08 = import_MSTmap(linkage_map_file_1)
individuals_01, genotypes_at_locus_01 = import_MSTmap(linkage_map_file_2)
individuals_09, genotypes_at_locus_09 = import_MSTmap(linkage_map_file_3)

my_pd_genos_08 = prep_data_pandas(individuals_08, genotypes_at_locus_08)
my_pd_genos_01 = prep_data_pandas(individuals_01, genotypes_at_locus_01)
my_pd_genos_09  = prep_data_pandas(individuals_09,  genotypes_at_locus_09)

all_my_data, loci_all = prepare_matrix(my_pd_genos_08, my_pd_genos_01, my_pd_genos_09)


fam_08, loci_08 = prepare_matrix(my_pd_genos_08)
fam_01, loci_01 = prepare_matrix(my_pd_genos_01)
fam_09, loci_09 = prepare_matrix(my_pd_genos_09)

# fam is now a tuple -> (genotype_matrix, loci)
fam_08 = prepare_matrix(my_pd_genos_08)
fam_01 = prepare_matrix(my_pd_genos_01)
fam_09 = prepare_matrix(my_pd_genos_09)

# rename markers here
def rename_loci_by_family(paralogs_file, fam_names, families):
    # check if each family listed in families is formatted as if returned from prepare matrix
    for afam in families:
        if not isinstance(afam, pd.core.frame.DataFrame):
            raise ValueError("families should be a pandas.core.frame.DataFrame")
        else:
            pass
    if len(fam_names) != len(families) :
        raise ValueError("names and families shoul have the same length")
    if not isinstance(fam_names, list ):
        raise ValueError("names should be a list")
    
    with open(blacklist_file) as x: 
        paralogs = [yy.strip() for yy in x.readlines()]
    # for each family, for each locus, if the locus is a paralog append family-specific text to locus name
    # genotypes are unchanged
    #new_familes = list()
    for idx, afam in enumerate(families):
        old_locus_names = afam.columns.values.tolist()
        new_locus_names = []
        for xx in old_locus_names:
            base_name = xx[:-3]
            if base_name in paralogs:
                print("{} is a paralog".format(base_name))
                new_name = "{}_{}_{}".format(base_name, fam_names[idx], xx[-2:])
            else: 
                new_name = base_name
            new_locus_names.append(new_name)
        afam.columns = new_locus_names
    return(families)
        
                
renamed_08, renamed_01, renamed_09 = rename_loci_by_family(paralogs_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_paralogs.txt", 
    fam_names = ['chum_08', 'chum_01', 'chum_09'], families = [my_pd_genos_08, my_pd_genos_01, my_pd_genos_09])

aa = pd.merge(left = renamed_08.transpose(), right = renamed_01.transpose(), how = 'outer', left_index =True, right_index =True)
bb = pd.merge(left = aa, right = renamed_09.transpose(), how = 'outer', left_index =True, right_index =True)  
bb.transpose()
all_my_data, loci_all = prepare_matrix(bb.transpose())


fams = [individuals_08, individuals_01, individuals_09]
LEPmap_filename = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap"
my_genotypes = bb.transpose()
my_genotypes = my_genotypes.replace(to_replace = [np.NaN, 0, 1, 2 ], value = ['0 0', '0 0', '1 1', '1 2'])

def write_LEPmap(families, family_names, loci, genotypes, output_filename):
    with open(output_filename, 'w') as OUTFILE:
        # TODO: write header
        header = "\t".join(["#family", 'name', 'sire', 'dam', 'sex', 'blank'] + loci) + "\n"
        OUTFILE.write(header)
        for fam_idx, fam in enumerate(families):
            fam_name = family_names[fam_idx]
            DAM_line = "\t".join([fam_name, fam_name + "_Dam", '0', '0', '2', '0'] + ['1 1' for xx in loci]) + "\n"
            SIRE_line = "\t".join([fam_name, fam_name + "_Sire", '0', '0', '1', '0'] + ['1 2' for xx in loci]) + "\n"
            OUTFILE.write(DAM_line)
            OUTFILE.write(SIRE_line)
            for ind in fam:
                ind_info = "\t".join([fam_name, ind, fam_name + "_Sire", fam_name + "_Dam", '0', '0'])
                ind_genotypes = genotypes.loc[ind]
                OUTFILE.write(ind_info + "\t" + "\t".join([str(xx) for xx in ind_genotypes]) + "\n")

write_LEPmap(families = fams, family_names = ["fam_08", "fam_01", "fam_09"], loci = loci_all, genotypes = my_genotypes, output_filename = LEPmap_filename)

print "java SeparateChromosomes data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
lodLimit = 10 sizeLimit = 3 \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes.log\
"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes \
lodLimit = 8 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod8_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod8_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod8_singles.chromosomes \
lodLimit = 7 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod7_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod7_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod7_singles.chromosomes \
lodLimit = 6 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod6_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod6_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod6_singles.chromosomes \
lodLimit = 5 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod5_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod5_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod5_singles.chromosomes \
lodLimit = 4 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod4_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod4_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod4_singles.chromosomes \
lodLimit = 3.5 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod3-5_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod3-5_singles.chromosomes.log"


all_fams = get_recombination_stats(all_my_data)
write_rec_stats(all_fams, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/all")

# Remove blacklisted markers (duplicates):
with open(blacklist_file) as x: 
    my_blacklist = [yy.strip() for yy in x.readlines()]
remove_by_blacklist(my_blacklist, genotypes_at_locus_08)
remove_by_blacklist(my_blacklist, genotypes_at_locus_01)
remove_by_blacklist(my_blacklist, genotypes_at_locus_09)

my_pd_genos_08 = prep_data_pandas(individuals_08, genotypes_at_locus_08)
my_pd_genos_01 = prep_data_pandas(individuals_01, genotypes_at_locus_01)
my_pd_genos_09  = prep_data_pandas(individuals_09,  genotypes_at_locus_09)

all_my_data, loci_all = prepare_matrix(my_pd_genos_08, my_pd_genos_01, my_pd_genos_09)

# write loci to file
write_loci(loci = loci_08, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_08")
write_loci(loci = loci_01, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_01")
write_loci(loci = loci_09, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_09")
write_loci(loci = loci_all, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/all")

#Each family alone
just_fam_08 = get_recombination_stats(fam_08)
just_fam_01 = get_recombination_stats(fam_01)
just_fam_09 = get_recombination_stats(fam_09)

write_rec_stats(just_fam_08, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_08")
write_rec_stats(just_fam_01, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_01")
write_rec_stats(just_fam_09, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats//fam_09")

# all families combined
all_fams = get_recombination_stats(all_my_data)
write_rec_stats(all_fams, "/olympus/WORK/WAPLES/chum_populations/rec_stats/all")

all_fams_minus_14 = get_recombination_stats(all_minus_14)
write_rec_stats(all_fams_minus_14, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all_minus_14")

# write LEPmap output
#loci_all # these are the 
fams = [individuals_08, individuals_01, individuals_09]
LEPmap_filename = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap"
my_genotypes = pd.concat(objs= [my_pd_genos_08, my_pd_genos_01,my_pd_genos_09])
my_genotypes = my_genotypes.replace(to_replace = [np.NaN, 0, 1, 2 ], value = ['0 0', '0 0', '1 1', '1 2'])

def write_LEPmap(families, family_names, loci, genotypes, output_filename):
    with open(output_filename, 'w') as OUTFILE:
        # TODO: write header
        header = "\t".join(["#family", 'name', 'sire', 'dam', 'sex', 'blank'] + loci) + "\n"
        OUTFILE.write(header)
        for fam_idx, fam in enumerate(families):
            fam_name = family_names[fam_idx]
            DAM_line = "\t".join([fam_name, fam_name + "_Dam", '0', '0', '2', '0'] + ['1 1' for xx in loci]) + "\n"
            SIRE_line = "\t".join([fam_name, fam_name + "_Sire", '0', '0', '1', '0'] + ['1 2' for xx in loci]) + "\n"
            OUTFILE.write(DAM_line)
            OUTFILE.write(SIRE_line)
            for ind in fam:
                ind_info = "\t".join([fam_name, ind, fam_name + "_Sire", fam_name + "_Dam", '0', '0'])
                ind_genotypes = genotypes.loc[ind]
                OUTFILE.write(ind_info + "\t" + "\t".join([str(xx) for xx in ind_genotypes]) + "\n")

write_LEPmap(families = fams, family_names = ["fam_08", "fam_01", "fam_09"], loci = loci_all, genotypes = my_genotypes, output_filename = LEPmap_filename)



# LEPmap commands
"java SeparateChromosomes data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.chromosomes"

"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.chromosomes \
lodLimit = 7 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes"

# chum_08


"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.chromosomes \
lodLimit = 7 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes"  

"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.chromosomes \
lodLimit = 7 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes"

for xx in [5,35]:
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
chromosome={} \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chr_{}.map \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chr_{}.map.log \
".format(xx, xx, xx)
) 

    
#estimate LOD limit

#ordermarkers
"java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
chromosome=1 \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_1.map \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_1.map.log \
" 

for xx in reversed(range(42)):
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
        chromosome={} \
        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map \
        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map.log\n".format(xx, xx, xx)
        )

for xx in [5,35]:
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
        chromosome={} \
        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map \
        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map.log\n".format(xx, xx, xx)
        )



# inspect
all_my_data.shape
len(genotypes_at_locus_14)








N = testy.R + testy.NR
R = testy.R
testa = e**-(2*(N/2. - R)**2/N)






#####################################
both_genos = pd.concat(objs = [my_pd_genos_14, my_pd_genos_10, my_pd_genos_A, my_pd_genos_B], join = 'outer')
# convert Nan to 0
both_genos.fillna(value = 0, inplace = True)
# conver to integer numpy array
both_genos_np = np.array(both_genos)
both_genos_np = both_genos_np.astype(int)
# transpose it
both_genos_np = np.transpose(both_genos_np)
#my_inds_14, my_genos_14 = import_MSTmap(linkage_map_file_1)
#my_inds_10, my_genos_10 = import_MSTmap(linkage_map_file_2)
#my_inds_A, my_genos_A = import_MSTmap(linkage_map_file_3)
#my_inds_B, my_genos_B = import_MSTmap(linkage_map_file_4)


######################
# Pandas treatment, under construction
# NOTICE, this is transposed from the numpy appraoch
# columns are loci, rows are individuals
my_pd_genos_14 = pd.DataFrame.from_dict(my_genos_14)
my_pd_genos_10 = pd.DataFrame.from_dict(my_genos_10)
my_pd_genos_A = pd.DataFrame.from_dict(my_genos_A)
my_pd_genos_B = pd.DataFrame.from_dict(my_genos_B)

# add an index
my_pd_genos_14.index = my_inds_14
my_pd_genos_10.index = my_inds_10
my_pd_genos_A.index = my_inds_A
my_pd_genos_B.index = my_inds_B




both_genos_np.shape


######################
# numpy treatment, WORKS! but no way to combine yet
#int_arr = np.array(my_genos.values())
# subset to just 100 loci
int_arr = 


