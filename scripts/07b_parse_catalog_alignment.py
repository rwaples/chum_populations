import pandas as pd
import numpy as np

# SAM file without header lines
SAM_filename = "C:/Users/IPSEG/Desktop/Waples/Dropbox/Projects/chum_populations/data/stacks_catalog/mapping/batch_42_self.sam"

# takes a basic SAM file and removes lines representing the self alignment (readName == refName)
with open(SAM_filename) as INFILE:
    with open(SAM_filename + "_self_removed", 'w') as OUTFILE:
        for line in INFILE:
            if line.startswith("@"):
                pass
            else:
                (readName, sumFlags, refName, refPos, mapQ, CIGAR, mate, matePos, fragmentLen, readSeq, readQ, alignScore), otherScores  = line.strip().split("\t")[:12], line.strip().split("\t")[12:] 
                if readName != refName:
                    OUTFILE.write(line)
                    
SAM_self_removed_filename = SAM_filename + "_self_removed"
whitelist_file = SAM_filename + "_whitelist"
blacklist_file = SAM_filename + "_blacklist"
mapped_loci_filename = "C:/Users/IPSEG/Desktop/Waples/Dropbox/Projects/chum_populations/data/stacks_catalog/mapping/mapped.stats.txt"

SAM = pd.read_csv(SAM_self_removed_filename, index_col = 0, sep = "\t", header = None)
SAM.index.names = ['catID']
SAM.columns = ["sumFlags", "refName", "refPos", "mapQ", "CIGAR", "mate", "matePos", "fragmentLen", "readSeq", "readQ", "alignScore", "XS", "XN", "XM", "XO", "XG", "NM", "MD", "YT"]
SAM.alignScore = [np.int(xx.split(":")[2]) for xx in SAM.alignScore]
SAM.head()
mapped_loci_dataframe = pd.read_csv(mapped_loci_filename, index_col = 'catID', sep = "\t", header = 0)
mapped_loci_dataframe.head()
mapped_loci = list()
for xx in mapped_loci_dataframe.index:
    try:
        mapped_loci.append(np.int(xx))
    except:
        pass
        
whitelist = set()

whitelist.update(mapped_loci)

blacklist = set()
with open(SAM_self_removed_filename) as INFILE:
    for line in INFILE:
        if line.startswith("@"):
            pass
        else:
            readName, sumFlags, refName, refPos, mapQ, CIGAR, mate, matePos, fragmentLen, readSeq, readQ, alignScore, XS, XN, XM, XO, XG, NM, MD, YT = line.strip().split("\t")
            readName = np.int(readName)
            refName = np.int(refName)
            if readName in whitelist and refName in whitelist:
                pass
            elif refName in whitelist and np.int(NM.split(":")[2]) <= 5:
                blacklist.add(np.int(readName))
            elif readName in whitelist and np.int(NM.split(":")[2]) <= 5:
                blacklist.add(np.int(refName))
            elif readName not in whitelist and refName not in whitelist and readName not in blacklist:
                whitelist.add(np.int(readName))
                blacklist.add(np.int(refName))

len(whitelist)
len(blacklist)
len(whitelist.intersection(blacklist))
    
with open(blacklist_file, 'w') as BLACKLIST:
   BLACKLIST.writelines([(str(x)+"\n") for x in blacklist])

# filter consensus FASTA with blacklist

with open(blacklist_file) as BLACKLIST:
    blacklist = [x.strip() for x in BLACKLIST.readlines()]
    
with open("C:/Users/IPSEG/Desktop/Waples/Dropbox/Projects/chum_populations/data/stacks_catalog/mapping/batch_42_consensus.fasta.txt") as INFILE:
    with open("C:/Users/IPSEG/Desktop/Waples/Dropbox/Projects/chum_populations/data/stacks_catalog/mapping/batch_42_CURATED.fasta.txt", 'w') as OUTFILE:
        for line in INFILE:
            if line.startswith(">"):
                if line.strip().replace(">", "") in blacklist:
                    print("excluding {}".format(line))
                else:
                    OUTFILE.write(line)
                    OUTFILE.write(next(INFILE))

# also added the chum salmon mtDNA genome (GenBank ID = AP010773 AB205144)
        


