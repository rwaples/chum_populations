from datetime import datetime
import os
import collections
import numpy as np
import itertools
import scipy.spatial.distance
#import operator
import pandas as pd
import sys

os.chdir("Y:/WORK/WAPLES/Stacks_mapping/Python")
from switch_allele_functions import (
    #parse_map_file_MST, 
    #convert_genotypes_to_int_array, 
    getR, getNR, 
    get_ml_R_frac, 
    get_LOD)


# Settings
linkage_map_file_1 = "Y:/WORK/MCKINNEY/MSTMap/MST_input/KUWHap14_mstmap_final.txt"
linkage_map_file_2 = "Y:/WORK/MCKINNEY/MSTMap/MST_input/Hap10_mstmap_final.txt"
linkage_map_file_3 = "Y:/WORK/MCKINNEY/MSTMap/MST_input/HapA_mstmap_final.txt"
linkage_map_file_4 = "Y:/WORK/MCKINNEY/MSTMap/MST_input/HapB_mstmap_final.txt"
blacklist_file = "Y:/WORK/MCKINNEY/MSTMap/blacklists/chinook_duplicates.txt"


#Functions
def import_MSTmap(filename):
    with open(filename) as INFILE:
        genotypes_at_locus = dict()
        for line in INFILE:
            if line.startswith('locus_name'): # start parsing
                individuals = line.strip().split()[1:]
                #print(individuals)
                for line in INFILE:
                    line_split = line.strip().split()
                    locus = line_split[0]
                    #print(locus)
                    genotypes = [1 if xx == 'a' else 2 if xx == 'b' else 0 for xx in line_split[1:]]
                    genotypes_at_locus[locus] = genotypes
                    #print(genotypes)
            else:
                pass
    return(individuals, genotypes_at_locus)
    
def remove_by_blacklist(blacklist, genotypes_at_locus):
    print("Starting length of genotypes: {}".format(len(genotypes_at_locus)))
    for locus in blacklist:        
        if locus in genotypes_at_locus.keys():
            del genotypes_at_locus[locus]
        elif locus + "_x1" in genotypes_at_locus.keys():
            del genotypes_at_locus[locus + "_x1"]
        elif locus + "_x2" in genotypes_at_locus.keys():
            del genotypes_at_locus[locus + "_x2"]
    print("Final length of genotypes: {}".format(len(genotypes_at_locus)))
    return(True)
    
def prep_data_pandas(individuals, genotypes_at_locus):
    # prepare pandas data.frame
    # columns are loci, rows are individuals
    my_pd_genos = pd.DataFrame.from_dict(genotypes_at_locus)
    # add an index
    my_pd_genos.index = individuals
    return(my_pd_genos)

def prepare_matrix(*args):
    for arg in args:
        if not isinstance(arg, pd.DataFrame):
            raise ValueError("need a pandas DataFrame")
        else:
            pass
    all_genos = pd.concat(objs = args, join = 'outer')
    # convert Nan to 0
    all_genos.fillna(value = 0, inplace = True)
    
    loci = [str(xx) for xx in all_genos.transpose().index]
    
    # conver to integer numpy array
    all_genos_np = np.array(all_genos)
    all_genos_np = all_genos_np.astype(int)
    # transpose it
    all_genos_np = np.transpose(all_genos_np)
    return(all_genos_np, loci)
    

# returns a redundant square matix
def get_matrix(data_in_array):
    data_in_matrix = scipy.spatial.distance.squareform(data_in_array)
    np.fill_diagonal(data_in_matrix, np.nan)
    return(data_in_matrix)

def get_recombination_stats(geno_array):
    int_arr = geno_array
    num_loci = int_arr.shape[0]
    num_pairs =  int((num_loci * (num_loci-1))/2)
    
    print('Starting, num_pairs = {}'.format(num_pairs))
    print(str(datetime.now()))
    time_start = datetime.now()
    sys.stdout.flush()
    
    pairs = itertools.combinations(int_arr, 2)
    R = np.fromiter(getR(pairs), dtype = np.int, count = num_pairs)
    time_R = datetime.now()
    print('Finished R')
    print(str(datetime.now()))
    sys.stdout.flush()
    
    pairs = itertools.combinations(int_arr, 2)
    NR = np.fromiter(getNR(pairs), dtype = np.int, count = num_pairs)
    time_NR = datetime.now()
    print('Finished NR')
    print(str(datetime.now()))
    sys.stdout.flush()
    
    ml_R_frac = get_ml_R_frac(R = R, NR = NR)
    time_RF = datetime.now()
    print('Finished RF')
    print(str(datetime.now()))
    
    sys.stdout.flush()
    Z = get_LOD(R = R, NR = NR, R_frac = ml_R_frac)
    time_Z = datetime.now()
    print('Finished Z')
    print(str(datetime.now()))
    sys.stdout.flush()
    
    N = R + NR
    MST = e**-(2*(N/2. - R)**2/N)
    print('Finished MST')
    time_MST = datetime.now()
    print(str(datetime.now()))
    sys.stdout.flush()    
    
    print("R took: {}".format(str(time_R - time_start)))
    print("NR took: {}".format(str(time_NR - time_R)))
    print("RF took: {}".format(str(time_RF - time_NR)))
    print("Z took: {}".format(str(time_Z - time_RF)))
    print("MST took: {}".format(str(time_MST - time_Z)))   
    
    Z_mat = get_matrix(Z)
    RF_mat = get_matrix(ml_R_frac)
    R_mat = get_matrix(R)
    NR_mat = get_matrix(NR)
    MST_mat = get_matrix(MST)
    
    Recombination_stats = collections.namedtuple('Recombination_stats', "R NR RF Z MST" )
    my_stats = Recombination_stats(R_mat, NR_mat, RF_mat, Z_mat, MST_mat)
    return(my_stats)    

def write_loci(loci, path):
    with open(os.path.join(path, 'loci.txt'), 'w') as OUTFILE:
        OUTFILE.write("\n".join(loci))
        
def write_rec_stats(stats, path):
    for stat in stats._fields:
        np.savetxt(X = getattr(stats, stat), fname = os.path.join(path, stat + ".tsv"), delimiter = "\t", fmt = '%1.4g')
            
# Import data from MST map input files:
individuals_14, genotypes_at_locus_14 = import_MSTmap(linkage_map_file_1)
individuals_10, genotypes_at_locus_10 = import_MSTmap(linkage_map_file_2)
individuals_A, genotypes_at_locus_A = import_MSTmap(linkage_map_file_3)
individuals_B, genotypes_at_locus_B = import_MSTmap(linkage_map_file_4)

my_pd_genos_14 = prep_data_pandas(individuals_14, genotypes_at_locus_14)
my_pd_genos_10 = prep_data_pandas(individuals_10, genotypes_at_locus_10)
my_pd_genos_A  = prep_data_pandas(individuals_A,  genotypes_at_locus_A)
my_pd_genos_B  = prep_data_pandas(individuals_B,  genotypes_at_locus_B)

fam_14, loci_14 = prepare_matrix(my_pd_genos_14)
fam_10, loci_10 = prepare_matrix(my_pd_genos_10)
fam_A, loci_A = prepare_matrix(my_pd_genos_A)
fam_B, loci_B = prepare_matrix(my_pd_genos_B)

# Remove blacklisted markers (duplicates):
with open(blacklist_file) as x: 
    my_blacklist = [yy.strip() for yy in x.readlines()]
remove_by_blacklist(my_blacklist, genotypes_at_locus_14)
remove_by_blacklist(my_blacklist, genotypes_at_locus_10)
remove_by_blacklist(my_blacklist, genotypes_at_locus_A)
remove_by_blacklist(my_blacklist, genotypes_at_locus_B)

my_pd_genos_14 = prep_data_pandas(individuals_14, genotypes_at_locus_14)
my_pd_genos_10 = prep_data_pandas(individuals_10, genotypes_at_locus_10)
my_pd_genos_A  = prep_data_pandas(individuals_A,  genotypes_at_locus_A)
my_pd_genos_B  = prep_data_pandas(individuals_B,  genotypes_at_locus_B)


all_my_data, loci_all = prepare_matrix(my_pd_genos_14, my_pd_genos_10, my_pd_genos_A, my_pd_genos_B)

all_minus_14, loci_all_minus_14 = prepare_matrix(my_pd_genos_10, my_pd_genos_A, my_pd_genos_B)



# write loci to file
write_loci(loci = loci_14, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_14")
write_loci(loci = loci_10, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_10")
write_loci(loci = loci_A, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_A")
write_loci(loci = loci_B, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_B")
write_loci(loci = loci_all, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all")
write_loci(loci = loci_all_minus_14, path = "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all_minus_14")


#Each family alone
just_fam_14 = get_recombination_stats(fam_14)
just_fam_10 = get_recombination_stats(fam_10)
just_fam_A = get_recombination_stats(fam_A)
just_fam_B = get_recombination_stats(fam_B)

write_rec_stats(just_fam_14, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_14")
write_rec_stats(just_fam_10, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_10")
write_rec_stats(just_fam_A, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_A")
write_rec_stats(just_fam_B, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/fam_B")

# all families combined
all_fams = get_recombination_stats(all_my_data)
write_rec_stats(all_fams, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all")

all_fams_minus_14 = get_recombination_stats(all_minus_14)
write_rec_stats(all_fams_minus_14, "Y:/WORK/MCKINNEY/MSTMap/recombination_stats/all_minus_14")


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










