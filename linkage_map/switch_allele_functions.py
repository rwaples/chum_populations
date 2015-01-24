import collections
import numpy
import itertools
import scipy.spatial.distance
import operator

def parse_map_file_rqtl(linkage_map_file):
    #this needs to return:
        # ini_map <- list of tuples(locus, lg, position)  *position unused  
        # loci_on_LG # count of loci on each LG (badly named)
    loci = list()
    lgs = list()
    loci_on_lg = collections.OrderedDict()
    positions = list() 
    with open (linkage_map_file, 'r') as INFILE:
        for line in INFILE:

            locus, lg, pos = line.strip().split()
            loci.append(locus)
            lgs.append(lg)
            if lg in loci_on_lg:
                loci_on_lg[lg] = loci_on_lg[lg]+1
            else:
                loci_on_lg[lg] = 1
            positions.append(pos)
    ini_map = zip(loci, lgs, positions)
    return(ini_map, loci_on_lg)
    
def parse_map_file_MST(linkage_map_file):
    line_count = 0
    loci = list()
    lgs = list()
    loci_on_lg = collections.OrderedDict()
    positions = list() 
    with open (linkage_map_file, 'r') as INFILE:
        current_LG = None
        for line in INFILE:
            line_count += 1
            if line.startswith(";"):
                pass #comment*ish lines
            elif line.startswith('group'):
                #set current_LG
                split_line = line.strip().split()
                lg_name = split_line[1]
                current_LG = lg_name[2:] # remove the preceeding text 'lg'
            elif not line.strip(): #test if empty string
                pass
            elif len(line.strip().split()) == 2:
                locus, pos = line.strip().split()
                loci.append(locus)
                lgs.append(current_LG)
                positions.append(pos)
                if current_LG in loci_on_lg:
                    loci_on_lg[current_LG] = loci_on_lg[current_LG]+1
                else:
                    loci_on_lg[current_LG] = 1
            else:
                raise ValueError("Unexpected Value at line {}".format(line_count))
    ini_map = zip(loci, lgs, positions)            
    return(ini_map, loci_on_lg)

def convert_genotypes_to_int_array(genotypes_of_locus, ini_map):
    mapped_genotypes = list()
    for locus_tuple in ini_map:
        locus_name = locus_tuple[0]
        mapped_genotypes.append(genotypes_of_locus[locus_name])
    geno_arr = numpy.array(mapped_genotypes)
    my_a = numpy.where(geno_arr == 'a')
    my_b = numpy.where(geno_arr == 'b')
    my_m = numpy.where(geno_arr == '-')
    int_arr = numpy.empty(shape = (geno_arr.shape), dtype = numpy.int)
    int_arr[my_a] = 1 # allele 1
    int_arr[my_b] = 2 # allele 2
    int_arr[my_m] = 0 # missing
    return(int_arr)


def getR(pairs):
    """generator, yields the numer of recombinant individuals for a given pair of loci"""
    for x, y in pairs:
        mult = x * y
        yield numpy.sum(mult == 2)
         
def getNR(pairs):
    """generator, yields the numer of NON-recombinant individuals for a given pair of loci"""
    for x, y in pairs:
        mult = x * y
        yield (numpy.sum(mult == 1) + numpy.sum(mult == 4))
    
def getM_iter(pairs):
    """generator, yields the numer of individuals missing calls for either of a given pair of loci"""
    for x, y in pairs:
        mult = x * y
        yield (sum(mult == 0))
    
def get_ml_R_frac(R, NR):
    """Returns the maximum likelyhood recombination fraction, as given as:
	R / (R + NR)"""
    if R.__class__ == numpy.ndarray:
	R = R.astype(numpy.float)	
    return(R / (R + NR))
    
def get_LOD(R, NR, R_frac):
    """Returns the LOD score for a set of the given values of [R, NR, R_frac].
    LOD scores calculated against the likelihood c =.5
    """
    Z = numpy.log10(
                (numpy.power((1-R_frac), NR) * numpy.power(R_frac, R)) / numpy.power(.5, (R + NR))
            )
    return(Z)
    
def get_threshold_recombinants_for_same_LGs(n, P):
    """
    Determine the threshold # of recombinations for a pair of loci to be (directly) placed on the same linkage group.
	n = total number of individuals 
	P = probability of false grouping of loci on same linkage group. Type I error. (unclear on interpretation)
	Same threshold used for all pairs of loci, not conditioned on missinging for the target pair.
    Used by MSTmap
    Taken from:
    Article Source:     
    Wu Y, Bhat PR, Close TJ, Lonardi S (2008) Efficient and Accurate Construction of Genetic Linkage Maps from the Minimum Spanning Tree of a Graph. PLoS Genet 4(10): e1000212. doi:10.1371/journal.pgen.1000212 Efficient and Accurate Construction of Genetic Linkage Maps from the Minimum Spanning Tree of a Graph
    """
    # If the number of NON-recombinants is below the number above then if the alleles were switch at one locus, 
    # they would be placed on the same linkage group.
    # NR < sig
    sig = n/2. - numpy.sqrt(n*numpy.log(P)/-2.)
    return(sig)
    
#Change to sqaure matrices
def get_rf_matrix(ml_R_frac):
    rf = scipy.spatial.distance.squareform(ml_R_frac)
    numpy.fill_diagonal(rf, numpy.nan)
    return(rf)

def get_lod_matrix(Z):
    lod = scipy.spatial.distance.squareform(Z)
    numpy.fill_diagonal(lod, numpy.nan)
    return(lod)
    
def get_NR_matrix(NR):
    NR_matrix = scipy.spatial.distance.squareform(NR)
    numpy.fill_diagonal(NR_matrix, numpy.nan)
    return(NR_matrix)
	
def convert_to_squareform(one_d_array):
    """Returns the supplied 1d array as a squareform distance matrix.
    Fills the diagonal with numpy.nan(s)"""
    two_d_array = scipy.spatial.distance.squareform(one_d_array)
    numpy.fill_diagonal(two_d_array, numpy.nan)
    return(two_d_array)


def get_index_of_LG(loci_on_lg):
    """
    Returns a dict of the (start, stop) positions used on 2d arrays to index the target LG.
    RELIES ON loci_on_lg BEING AN OrderedDict"""
    index_of_lg = collections.OrderedDict()
    index = 0
    for lg, num_loci in loci_on_lg.items():
        start = index
        stop = start + num_loci
        index = stop
        index_of_lg[lg] = start, stop
    return(index_of_lg)
    
def find_LGs_with_multiple_loci(index_of_lg, loci_on_lg):
    ordered_lgs = index_of_lg.keys()
    lgs_longer_than_1 = list()
    for xx in ordered_lgs:
        if loci_on_lg[xx] > 1:
            lgs_longer_than_1.append(xx)
    return(lgs_longer_than_1)

def get_mean(pairs, array, index_of_lg):
    """generator, yields the mean value of the array for each pair of LGs in paris"""
    for x, y in pairs:
        lg_a_start = index_of_lg[x][0]
        lg_a_stop = index_of_lg[x][1]
        lg_b_start = index_of_lg[y][0]
        lg_b_stop = index_of_lg[y][1]
        mean = numpy.mean(array[lg_a_start:lg_a_stop, lg_b_start:lg_b_stop].flatten())
        yield (mean)
        
def get_sum(pairs, array, index_of_lg):
    """generator, yields the sumof the array for each pair of LGs in paris"""
    for x, y in pairs:
        lg_a_start = index_of_lg[x][0]
        lg_a_stop = index_of_lg[x][1]
        lg_b_start = index_of_lg[y][0]
        lg_b_stop = index_of_lg[y][1]
        my_sum = numpy.sum(array[lg_a_start:lg_a_stop, lg_b_start:lg_b_stop].flatten())
        yield (my_sum)        
        
def get_count_NR_below_thresh(pairs, array, index_of_lg, threshold):
    """generator, yields the number of locus pairs that fall strictly below the theshold,
	for each pair of LGs in paris"""
    for x, y in pairs:
        lg_a_start = index_of_lg[x][0]
        lg_a_stop = index_of_lg[x][1]
        lg_b_start = index_of_lg[y][0]
        lg_b_stop = index_of_lg[y][1]
        target = array[lg_a_start:lg_a_stop, lg_b_start:lg_b_stop].flatten()
        my_sum = numpy.sum(target < threshold)
        yield (my_sum)

def get_LG_pairwise_count_NR_threshold(lgs_longer_than_1, NR_matrix, index_of_lg, threshold):
    pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2) 
    NR_under_threshold = numpy.fromiter(get_count_NR_below_thresh(pairs_of_lgs, NR_matrix, index_of_lg, threshold), dtype = numpy.int)
    return(NR_under_threshold)

def get_LG_pairwise_mean_rf(lgs_longer_than_1, rf, index_of_lg):
    pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2) 
    mean_rf = numpy.fromiter(get_mean(pairs_of_lgs, rf, index_of_lg), dtype = numpy.float)
    return(mean_rf)
    
def get_LG_pairwise_mean_lod(lgs_longer_than_1, lod, index_of_lg):
    pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2) 
    mean_lod = numpy.fromiter(get_mean(pairs_of_lgs, lod, index_of_lg), dtype = numpy.float)
    return(mean_lod)
    
def get_LG_pairwise_sum_lod(lgs_longer_than_1, lod, index_of_lg):
    pairs_of_lgs = itertools.combinations(lgs_longer_than_1, 2) 
    sum_lod = numpy.fromiter(get_sum(pairs_of_lgs, lod, index_of_lg), dtype = numpy.float)
    return(sum_lod)
   

#change to square_matrix, only containing LG with multiple loci
# can this be combined with the other similar function?
def get_square_form(one_d_array, lgs_longer_than_1):    
    sq_array = scipy.spatial.distance.squareform(one_d_array)
    upper_tri = numpy.triu_indices(len(lgs_longer_than_1))
    sq_array[upper_tri] = 0
    return(sq_array)


def find_LG_pairs_with_sum_lod_above(sq_sum_lod, threshold):
    where_pairs = numpy.where(sq_sum_lod > threshold)
    return(where_pairs)



#Select LG that appears in most switched pairs, if tied pick one
# Add that LG to switch list
# remove all pairs contianing that LG from the initial list
# repeat until no pairs remain, return switch list

def find_LGs_to_switch(lgs_longer_than_1, where_pairs):
    matched_pairs = zip(*where_pairs)
    LGs_to_switch = list()
    while len(matched_pairs) > 0:
        count_switched = collections.defaultdict(int)
        for ii, jj in matched_pairs:
            count_switched[ii] += 1
            count_switched[jj] += 1
        sorted_switched = sorted(count_switched.items(),key=operator.itemgetter(1), reverse = False)   
        most_common_LG_index = sorted_switched.pop()[0]
        LG_name = lgs_longer_than_1[most_common_LG_index]
        print("most common LG: {}".format(LG_name))
        LGs_to_switch.append(LG_name)
        to_remove = list()
        for index in range(len(matched_pairs)):
            xx, yy = matched_pairs[index]
            if most_common_LG_index == xx or  most_common_LG_index == yy:
                to_remove.append((xx, yy))
        for ii in to_remove:    
            matched_pairs.remove(ii)
            print("Removed: {}".format(ii))
    #translate back to LG names, and not positions as in where_pairs
    return(LGs_to_switch)

def find_loci_to_switch(ini_map, LGs_to_switch):
    list_of_loci_on_LG = collections.defaultdict(list)                                                            
    for locus, LG, pos in ini_map:
        #convert to str
        LG = str(LG)
        list_of_loci_on_LG[LG].append(locus)
    loci_to_switch = list()
    for LG in LGs_to_switch:
        loci_to_switch += list_of_loci_on_LG[LG]
    print("Returned {} loci on {} LGs".format(len(loci_to_switch), len(LGs_to_switch)))
    return(loci_to_switch)
    
####################3    
    
def calculate_switch_stats(mappable, linkage_map_file, linkage_map_format, MST_grouping_threshold):
    genotypes_of_locus = mappable
    if linkage_map_format.lower() == 'mst':
        ini_map, loci_on_lg = parse_map_file_MST(linkage_map_file)
    elif linkage_map_format.lower() == 'rqtl':   
        ini_map, loci_on_lg = parse_map_file_rqtl(linkage_map_file)
    else:
        raise ValueError("unknown linkage_map_format")
    
    int_arr = convert_genotypes_to_int_array(genotypes_of_locus, ini_map)
    num_loci = int_arr.shape[0]
    num_pairs =  int((num_loci * (num_loci-1))/2)
    pairs = itertools.combinations(int_arr, 2)
    R = numpy.fromiter(getR(pairs), dtype = numpy.float64, count = num_pairs)
    pairs = itertools.combinations(int_arr, 2)
    NR = numpy.fromiter(getNR(pairs), dtype = numpy.float64, count = num_pairs)
    ml_R_frac = get_ml_R_frac(R = R, NR = NR)
    Z = get_LOD(R = R, NR = NR, R_frac = ml_R_frac)
    NR_matrix = get_NR_matrix(NR)
    #rf = get_rf_matrix(ml_R_frac)
    lod = get_lod_matrix(Z)
    index_of_lg = get_index_of_LG(loci_on_lg)
    lgs_longer_than_1 = find_LGs_with_multiple_loci(index_of_lg, loci_on_lg)
    #mean_rf = get_LG_pairwise_mean_rf(lgs_longer_than_1, rf, index_of_lg)
    #mean_lod = get_LG_pairwise_mean_lod(lgs_longer_than_1,lod, index_of_lg)
    sum_lod = get_LG_pairwise_sum_lod(lgs_longer_than_1,lod, index_of_lg)
    sq_sum_lod = get_square_form(sum_lod, lgs_longer_than_1)
    n = len(mappable.items()[0][1]) #number of individuals
    NR_threshold = get_threshold_recombinants_for_same_LGs(n, MST_grouping_threshold)
    NR_under_threshold = get_LG_pairwise_count_NR_threshold(lgs_longer_than_1, NR_matrix, index_of_lg, threshold = NR_threshold)
    sq_NR_matrix = get_square_form(NR_under_threshold, lgs_longer_than_1)
    return(ini_map, sq_sum_lod, sq_NR_matrix, R, NR, lgs_longer_than_1)
    
def switch_with_threshold(ini_map, sq_sum_lod, sum_lod_threshold, lgs_longer_than_1):
    where_pairs = find_LG_pairs_with_sum_lod_above(sq_sum_lod, sum_lod_threshold)
    LGs_to_switch = find_LGs_to_switch(lgs_longer_than_1, where_pairs, sq_sum_lod)
    loci_to_switch = find_loci_to_switch(ini_map, LGs_to_switch)
    return(LGs_to_switch, loci_to_switch)

#Actually do the allele switching        
def switch_alleles(loci_to_switch, mappable):
    for jj in loci_to_switch:
        if jj in mappable:
            original_genotypes = mappable[jj]
            switched_genotpyes = ['a' if xx == 'b' else 'b' if xx == 'a' else xx for xx in original_genotypes]
            mappable[jj] = switched_genotpyes
    return(mappable)
    
def update_switched_alleles_file(switched_alleles_file, loci_to_switch):
    set_loci_to_switch = set(loci_to_switch)
    pre_switched = set()
    try:
        with open(switched_alleles_file, 'r'): pass
    except IOError:
        try:
            with open(switched_alleles_file, 'w'): pass
        except IOError:
            print ("Bad file name")
    with open(switched_alleles_file, 'r') as infile:
        for line in infile:
            locus = line.strip().split()
            pre_switched.update(set([locus]))
    in_either = set_loci_to_switch.union(pre_switched)
    in_both = set_loci_to_switch.intersection(pre_switched)
    #write loci present in in_either but not in in_both
    loci_to_write = in_either.difference(in_both)
    print("Writing {} loci to switched alleles file".format(len(loci_to_write)))
    with open(switched_alleles_file, 'w') as outfile:
        for locus in loci_to_write:
            outfile.write("{}\n".format(locus))