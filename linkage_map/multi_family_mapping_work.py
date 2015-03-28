from multi_family_mapping_functions import *
import pandas as pd
import collections

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

#all_my_data, loci_all = prepare_matrix(my_pd_genos_08, my_pd_genos_01, my_pd_genos_09)


#fam_08, loci_08 = prepare_matrix(my_pd_genos_08)
#fam_01, loci_01 = prepare_matrix(my_pd_genos_01)
#fam_09, loci_09 = prepare_matrix(my_pd_genos_09)

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
    
    with open(paralogs_file) as INFILE: 
        paralogs = [yy.strip() for yy in INFILE.readlines()]
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


###### this should be moved inside function
renamed_08t = renamed_08.transpose()
renamed_01t = renamed_01.transpose()
renamed_09t = renamed_09.transpose()
aa = pd.merge(left = renamed_08t, right = renamed_01t, how = 'outer', left_index =True, right_index =True)
bb = pd.merge(left = aa, right = renamed_09t, how = 'outer', left_index =True, right_index =True) 
###### end move into function

all_my_data, loci_all = prepare_matrix(bb.transpose())

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


fams = [individuals_08, individuals_01, individuals_09]
LEPmap_filename = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap"
my_genotypes = bb.transpose()
my_genotypes = my_genotypes.replace(to_replace = [np.NaN, 0, 1, 2 ], value = ['0 0', '0 0', '1 1', '1 2'])
write_LEPmap(families = fams, family_names = ["fam_08", "fam_01", "fam_09"], loci = loci_all, genotypes = my_genotypes, output_filename = LEPmap_filename)

print "java SeparateChromosomes data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lepmap \
lodLimit = 10 sizeLimit = 3 \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes.log"

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

# this takes a long time!
#all_fams = get_recombination_stats(all_my_data)
#write_rec_stats(all_fams, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/all")
# end long time!


#look for LG assignments of duplicated loci that are congruent across families
def find_duplicate_names(paralogs_file, loci, LG_file, output_file):
    with open(blacklist_file) as INFILE: 
        paralogs = [yy.strip() for yy in INFILE.readlines()]
    #get list of the locus names
        #  given by loci
    #get list of LG assignments
    with open(LG_file) as INFILE:
        #skip first line
        next(INFILE)
        LG_assignments = [int(line.strip()) for line in INFILE]
    LG_of_locus = dict(zip(loci, LG_assignments))
    #print(paralogs)
    with open(output_file, 'w') as OUTFILE:
        OUTFILE.write("{}\t{}\t{}\t{}\n".format('catalog_name', 'family', 'copy', 'LG'))
        for locus in loci:
            catalog_name = locus.split("_")[0]
            if catalog_name in paralogs:
                family = locus.split("_")[2]
                copy = locus[-2:]
                OUTFILE.write("{}\t{}\t{}\t{}\n".format(catalog_name, family, copy, LG_of_locus[locus]))
            #list of all loci sharing base name
                sharing = [loc for loc in loci if locus.split("_")[0] == loc.split("_")[0]]
                sharing.remove(locus)
                agree_on_LG = sum([1 for loc in sharing if LG_of_locus[loc] == LG_of_locus[locus] ])
                disagree_on_LG = sum([1 for loc in sharing if LG_of_locus[loc] != LG_of_locus[locus] ])
                #print(locus, agree_on_LG, disagree_on_LG)

find_duplicate_names(paralogs_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_paralogs.txt", 
    loci = loci_all, LG_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.lod5_singles.chromosomes", 
    output_file = '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/LG_congruence.tsv')


def write_name_table(paralogs_file, loci, LG_file, output_file):
    with open(blacklist_file) as INFILE: 
        paralogs = [yy.strip() for yy in INFILE.readlines()]
    with open(LG_file) as INFILE:
        #skip first line
        next(INFILE)
        LG_assignments = [int(line.strip()) for line in INFILE]
    LG_of_locus = defaultdict(int, zip(loci, LG_assignments))
    with open(output_file, 'w') as OUTFILE:
        OUTFILE.write("{}\t{}\t{}\t{}\n".format('catalog_name', 'LG_08_x1', 'LG_08_x2', 'LG_01_x1', 'LG_01_x2', 'LG_09_x1', 'LG_09_x2'))
        for locus in loci:
            catalog_name = locus.split("_")[0]
            if catalog_name in paralogs:
                family = locus.split("_")[2]
                copy = locus[-2:]
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\n".format(catalog_name, LG_of_locus[locus]))
            #list of all loci sharing base name
                sharing = [loc for loc in loci if locus.split("_")[0] == loc.split("_")[0]]
                sharing.remove(locus)
                agree_on_LG = sum([1 for loc in sharing if LG_of_locus[loc] == LG_of_locus[locus] ])
                disagree_on_LG = sum([1 for loc in sharing if LG_of_locus[loc] != LG_of_locus[locus] ])
                #print(locus, agree_on_LG, disagree_on_LG)

find_duplicate_names(paralogs_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/chum_paralogs.txt", 
    loci = loci_all, LG_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/all_loci.inital_chromosomes", 
    output_file = '/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/name_table.tsv')


def collapse_names(LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2):
    # take them in order given here
    possible_results = ['A','B','C','D','E','F']
    #possible_results.reverse()
    result = []
    
    # check for segmental
    # if segmental do not try to resolve
    for xx in ((LG_08_x1, LG_08_x2), (LG_01_x1, LG_01_x2), (LG_09_x1, LG_09_x2)):
        x1, x2 = xx
        if x1 == x2 and x1 != 0: # segmental
            result = possible_results

    else:
        mapping = dict()
        for assign in (LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2):
            if assign in mapping:
                result.append(mapping[assign])
            else:
                if assign == 0:
                     result.append(possible_results.pop(0))
                else:
                    mapping[assign] = possible_results.pop(0)
                    result.append(mapping[assign])
    return(result)  
        #diffs = list(set(LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2))
        

    #if LG_08_x1 != 0:
    #    LG_A = LG_08_x1
    #    result.append("A")
    #else:
    #    result.append("NA")
    #    LG_A = None
    #if LG_08_x2 != 0:
    #    if result == ["NA"]:
    #        result.append("A")
    #        LG_A = LG_08_x2
    #        LG_B = None
    #    else:
    #        result.append("B")
    #        LG_B = LG_08_x2
    #        if LG_08_x1 == LG_08_x2:
    #            segmental == True
    #            result.extend(('NA', 'NA', 'NA', 'NA'))
    #            return(result)
    #else:
    #    result.append("NA")
    #    LG_B = None
    #    
    ## determine matching in other families
    #for xx in [LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2]:
    #    if xx == LG_A:
    #        result.append('A')
    #    elif xx == LG_B:
    #        result.append('B')
    #    elif xx != 0:
    #        if LG_A == None:
    #            LG_A = xx
    #            result.append('A')
    #        elif LG_B == None:
    #            LG_B = xx
    #            result.append('B')            
    #        else:
    #            result.append('C')
    #    else:
    #        result.append('NA')
    #return(result)

#examples
print collapse_names(1,2,2,1,2,1)
print collapse_names(0,1,2,0,2,0)
print collapse_names(1,0,2,0,2,0)
print collapse_names(0,0,1,0,2,0)
print collapse_names(0,12,0,0,2,2)

def parse_LG_congruence_line(line):
    catalog_name, family, copy, LG = line.strip().split("\t")
    return(catalog_name, family, copy, LG)


def write_rename_table(LG_congruence_file, out_file):
    famLG_of_locus = collections.defaultdict(dict)
    with open(out_file, 'w') as OUTFILE:
        OUTFILE.write("{}\t{}\t{}\n".format("old_name", "new_name", "LG"))
        with open(LG_congruence_file) as PARALOG_LGS:
            #skip first line
            next(PARALOG_LGS)
            for line in PARALOG_LGS:
                catalog_name, family, copy, LG = parse_LG_congruence_line(line)
                famLG_of_locus[catalog_name][family, copy] = int(LG)
            for cn, famLG in famLG_of_locus.items():
                LG_08_x1 = famLG.get(('08', 'x1'), 0)
                LG_08_x2 = famLG.get(('08', 'x2'), 0)
                LG_01_x1 = famLG.get(('01', 'x1'), 0)
                LG_01_x2 = famLG.get(('01', 'x2'), 0)
                LG_09_x1 = famLG.get(('09', 'x1'), 0)
                LG_09_x2 = famLG.get(('09', 'x2'), 0)
                #print cn, [LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2], collapse_names(LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2)
                fam_names = ["{}_chum_08_x1", "{}_chum_08_x2", "{}_chum_01_x1", "{}_chum_01_x2", "{}_chum_09_x1", "{}_chum_09_x2"]
                collaped_names = collapse_names(LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2)
                LGs = [LG_08_x1, LG_08_x2, LG_01_x1, LG_01_x2, LG_09_x1, LG_09_x2]
                family_specific_names = [xx.format(cn) for xx in fam_names]
                for cnt, fsn in enumerate(family_specific_names):
                    OUTFILE.write("{}\t{}\t{}\n".format(fsn, cn+"_{}".format(collaped_names[cnt]), LGs[cnt]))
        
write_rename_table(LG_congruence_file =  "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/LG_congruence.tsv", 
    out_file = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/rename_table.tsv")   

#get rename dict
rename_table = pd.read_table("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/rename_table.tsv")
rename_dict = dict(zip(rename_table.old_name,rename_table.new_name))
    
renamed_08t = renamed_08.rename(columns=rename_dict).transpose()
renamed_01t = renamed_01.rename(columns=rename_dict).transpose()
renamed_09t = renamed_09.rename(columns=rename_dict).transpose()
aa = pd.merge(left = renamed_08t, right = renamed_01t, how = 'outer', left_index =True, right_index =True)
bb = pd.merge(left = aa, right = renamed_09t, how = 'outer', left_index =True, right_index =True)

fams = [individuals_08, individuals_01, individuals_09]
LEPmap_filename = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap"
my_genotypes = bb.transpose()
my_genotypes = my_genotypes.replace(to_replace = [np.NaN, 0, 1, 2 ], value = ['0 0', '0 0', '1 1', '1 2'])
write_LEPmap(families = fams, family_names = ["fam_08", "fam_01", "fam_09"], loci = my_genotypes.columns.values.tolist(),
    genotypes = my_genotypes, output_filename = LEPmap_filename)

# form linkage groups
print "java SeparateChromosomes data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
lodLimit = 10 sizeLimit = 20 \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.inital.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.inital.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.inital.chromosomes \
lodLimit = 8 lodDifference=3 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod8_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod8_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod8_singles.chromosomes \
lodLimit = 7 lodDifference=3 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod7_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod7_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod7_singles.chromosomes \
lodLimit = 6 lodDifference=3 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod6_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod6_singles.chromosomes.log"

print "java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod6_singles.chromosomes \
lodLimit = 5 lodDifference=3 \
data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes \
2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes.log"


with open("/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes") as LG_file:
    next(LG_file)
    print zip(my_genotypes.columns.values.tolist(), [xx.strip() for xx in LG_file])



for xx in reversed(range(38,42)):
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes \
        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
        alpha=1 maxDistance=30 \
        chromosome={} \
        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map \
        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map.log\n".format(xx, xx, xx)
        )

for xx in (1, 7):
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes \
        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
        alpha=1 maxDistance=30 \
        chromosome={} \
        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map \
        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map.log\n".format(xx, xx, xx)
        )

for xx in (1, 7, 25, 26, 31, 33):
    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lod5_singles.chromosomes \
        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed_loci.lepmap \
        alpha=.5 maxDistance=20 \
        chromosome={} \
        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map \
        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/chr_{}.map.log\n".format(xx, xx, xx)
        )


#
######################################################################################################3
###END WORK
#
#
#
## different ways or reasons to rename
## renaming only occuring *within* catalog entries, and seeks to standardize naming across families
## non paralogs
#    # should already share names across families
#
## paralogs
#    # rename relative to a specified target family 
#    # if same catalog locus assigned to same LG in two different families, collapse names
#    # if same catalog locus assigned to same LG in the same families, do not collapse names, possible segmental duplicate
#    
#    # if same catalog locus assigned to different LG in two different families, check if the LGs are homeologs
#        # if homeologs - do not callpase names incormporating LG
#        # if not homeologs - still rename??
#        
## what to do when a locus is unmapped in the small families?
#
## renaming procedure
## 1 - identify homeologs based on family 8
## 2 - 
#
#    
#
## mapped two different families, same LG assignmentin small family 
#    
#    
#        
#        
#    
#'10001_y' in my_genotypes.columns.values.tolist()
#
#
## Remove blacklisted markers (duplicates):
#with open(blacklist_file) as x: 
#    my_blacklist = [yy.strip() for yy in x.readlines()]
#remove_by_blacklist(my_blacklist, genotypes_at_locus_08)
#remove_by_blacklist(my_blacklist, genotypes_at_locus_01)
#remove_by_blacklist(my_blacklist, genotypes_at_locus_09)
#
#my_pd_genos_08 = prep_data_pandas(individuals_08, genotypes_at_locus_08)
#my_pd_genos_01 = prep_data_pandas(individuals_01, genotypes_at_locus_01)
#my_pd_genos_09  = prep_data_pandas(individuals_09,  genotypes_at_locus_09)
#
#all_my_data, loci_all = prepare_matrix(my_pd_genos_08, my_pd_genos_01, my_pd_genos_09)
#
## write loci to file
#write_loci(loci = loci_08, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_08")
#write_loci(loci = loci_01, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_01")
#write_loci(loci = loci_09, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_09")
#write_loci(loci = loci_all, path = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/all")
#
##Each family alone
#just_fam_08 = get_recombination_stats(fam_08)
#just_fam_01 = get_recombination_stats(fam_01)
#just_fam_09 = get_recombination_stats(fam_09)
#
#write_rec_stats(just_fam_08, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_08")
#write_rec_stats(just_fam_01, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats/fam_01")
#write_rec_stats(just_fam_09, "/home/ipseg/Desktop/waples/chum_populations/linkage_map/rec_stats//fam_09")
#
## all families combined
#all_fams = get_recombination_stats(all_my_data)
#write_rec_stats(all_fams, "/olympus/WORK/WAPLES/chum_populations/rec_stats/all")
#
#all_fams_minus_14 = get_recombination_stats(all_minus_14)
#write_rec_stats(all_fams_minus_14, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all_minus_14")
#
## write LEPmap output
##loci_all # these are the 
#fams = [individuals_08, individuals_01, individuals_09]
#LEPmap_filename = "/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap"
#my_genotypes = pd.concat(objs= [my_pd_genos_08, my_pd_genos_01,my_pd_genos_09])
#my_genotypes = my_genotypes.replace(to_replace = [np.NaN, 0, 1, 2 ], value = ['0 0', '0 0', '1 1', '1 2'])
#
#def write_LEPmap(families, family_names, loci, genotypes, output_filename):
#    with open(output_filename, 'w') as OUTFILE:
#        # TODO: write header
#        header = "\t".join(["#family", 'name', 'sire', 'dam', 'sex', 'blank'] + loci) + "\n"
#        OUTFILE.write(header)
#        for fam_idx, fam in enumerate(families):
#            fam_name = family_names[fam_idx]
#            DAM_line = "\t".join([fam_name, fam_name + "_Dam", '0', '0', '2', '0'] + ['1 1' for xx in loci]) + "\n"
#            SIRE_line = "\t".join([fam_name, fam_name + "_Sire", '0', '0', '1', '0'] + ['1 2' for xx in loci]) + "\n"
#            OUTFILE.write(DAM_line)
#            OUTFILE.write(SIRE_line)
#            for ind in fam:
#                ind_info = "\t".join([fam_name, ind, fam_name + "_Sire", fam_name + "_Dam", '0', '0'])
#                ind_genotypes = genotypes.loc[ind]
#                OUTFILE.write(ind_info + "\t" + "\t".join([str(xx) for xx in ind_genotypes]) + "\n")
#
#write_LEPmap(families = fams, family_names = ["fam_08", "fam_01", "fam_09"], loci = loci_all, genotypes = my_genotypes, output_filename = LEPmap_filename)
#
#
#
## LEPmap commands
#"java SeparateChromosomes data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.chromosomes"
#
#"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.chromosomes \
#lodLimit = 7 \
#data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes"
#
## chum_08
#
#
#"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.chromosomes \
#lodLimit = 7 \
#data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes"  
#
#"java JoinSingles /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.chromosomes \
#lodLimit = 7 \
#data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes"
#
#for xx in [5,35]:
#    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap.final.chromosomes \
#data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chum_08.lepmap \
#chromosome={} \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chr_{}.map \
#2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chum_08/chr_{}.map.log \
#".format(xx, xx, xx)
#) 
#
#    
##estimate LOD limit
#
##ordermarkers
#"java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
#data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#chromosome=1 \
#> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_1.map \
#2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_1.map.log \
#" 
#
#for xx in reversed(range(42)):
#    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
#        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#        chromosome={} \
#        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map \
#        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map.log\n".format(xx, xx, xx)
#        )
#
#
#
#for xx in reversed(range(42)):
#    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
#        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#        chromosome={} \
#        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map \
#        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map.log\n".format(xx, xx, xx)
#        )
#
#for xx in [5,35]:
#    print("java OrderMarkers /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap.final.chromosomes \
#        data=/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/ini.lepmap \
#        chromosome={} \
#        > /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map \
#        2> /home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/chr_{}.map.log\n".format(xx, xx, xx)
#        )
#
#
#
## inspect
#all_my_data.shape
#len(genotypes_at_locus_14)
#
#
#
#
#
#
#
#
#N = testy.R + testy.NR
#R = testy.R
#testa = e**-(2*(N/2. - R)**2/N)
#
#
#
#
#
#
######################################
#both_genos = pd.concat(objs = [my_pd_genos_14, my_pd_genos_10, my_pd_genos_A, my_pd_genos_B], join = 'outer')
## convert Nan to 0
#both_genos.fillna(value = 0, inplace = True)
## conver to integer numpy array
#both_genos_np = np.array(both_genos)
#both_genos_np = both_genos_np.astype(int)
## transpose it
#both_genos_np = np.transpose(both_genos_np)
##my_inds_14, my_genos_14 = import_MSTmap(linkage_map_file_1)
##my_inds_10, my_genos_10 = import_MSTmap(linkage_map_file_2)
##my_inds_A, my_genos_A = import_MSTmap(linkage_map_file_3)
##my_inds_B, my_genos_B = import_MSTmap(linkage_map_file_4)
#
#
#######################
## Pandas treatment, under construction
## NOTICE, this is transposed from the numpy appraoch
## columns are loci, rows are individuals
#my_pd_genos_14 = pd.DataFrame.from_dict(my_genos_14)
#my_pd_genos_10 = pd.DataFrame.from_dict(my_genos_10)
#my_pd_genos_A = pd.DataFrame.from_dict(my_genos_A)
#my_pd_genos_B = pd.DataFrame.from_dict(my_genos_B)
#
## add an index
#my_pd_genos_14.index = my_inds_14
#my_pd_genos_10.index = my_inds_10
#my_pd_genos_A.index = my_inds_A
#my_pd_genos_B.index = my_inds_B
#
#
#
#
#both_genos_np.shape
#
#
#######################
## numpy treatment, WORKS! but no way to combine yet
##int_arr = np.array(my_genos.values())
## subset to just 100 loci
#int_arr = 
#

