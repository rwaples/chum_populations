from datetime import datetime
import os
import collections
import numpy as np
import itertools
import scipy.spatial.distance
import pandas as pd
import sys
from switch_allele_functions import (
    #parse_map_file_MST, 
    #convert_genotypes_to_int_array, 
    getR, getNR, 
    get_ml_R_frac, 
    get_LOD)


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
    MST = np.e**-(2*(N/2. - R)**2/N)
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

            
  







