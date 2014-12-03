import os
import os.path
import gzip
import glob
import struct


target_dir = "/media/Shared/Data/chum/populations/fastQC"


fastQC_list = glob.glob(os.path.join(target_dir, '*.fq_fastqc'))

def parse_name_fqgz(ind_filename):
    basename = os.path.basename(ind_filename)
    root, ext = os.path.splitext(basename)
    sample_name, fq = os.path.splitext(root)
    silli, ind_num = sample_name.split('_')
    return(basename, sample_name, silli, ind_num)
    
with open("/home/ipseg/Desktop/waples/chum_populations/results/ind_seq_stats.tsv", 'w') as OUTFILE:
    OUTFILE.write("\t".join(['filename', 'name', 'pop', 'ind', 'encoding', 'nreads', 'lenreads']))
    OUTFILE.write("\n")
    for fastQC_results_dir in fastQC_list:
        with open(os.path.join(fastQC_results_dir, "fastqc_data.txt")) as INFILE:
            for line in INFILE:
                line = line.strip()
                if line.startswith('Filename'):
                    sample_filename = line.split("\t")[1]
                elif line.startswith('Encoding'):
                    encoding = line.split("\t")[1]
                elif line.startswith('Total Sequences'):
                    nreads = line.split("\t")[1]
                elif line.startswith('Sequence length'):
                    lenreads = line.split("\t")[1]
            OUTFILE.write("\t".join(parse_name_fqgz(sample_filename)))
            OUTFILE.write("\t" + encoding)
            OUTFILE.write("\t" + nreads)
            OUTFILE.write("\t" + lenreads + "\n")          
            

            

### OLD WAY ###
fqgz_list = glob.glob(os.path.join("/media/Shared/Data/chum/populations/cleanSeqs", '*.fq.gz'))
#fqgz_list += glob.glob(os.path.join("/media/Shared/Data/chum/populations/cleanSeqs/PE_Hoodsport", '*.fq.gz'))



  
def gzipFileSize(filename):
    "return UNCOMPRESSED filesize of a gzipped file"
    fo = open(filename, 'rb')
    fo.seek(-4, 2)
    r = fo.read()
    fo.close()
    return struct.unpack('<I', r)[0]
    
def nlines_gzip(gzip_filename):
    with gzip.open(gzip_filename) as INFILE:
        nlines = sum(1 for ln in INFILE)
        return(nlines)
    
def lenreads_gzip(gzip_filename):
    with gzip.open(gzip_filename) as INFILE:
        first_line = next(INFILE)
        second_line = next(INFILE)
        return(len(second_line.strip())) 
    

with open("/home/ipseg/Desktop/waples/chum_populations/results/ind_seq_stats.tsv", 'w') as OUTFILE:
    OUTFILE.write("\t".join(['filename', 'name', 'pop', 'ind', 'size_bytes', 'nreads', 'lenreads']))
    OUTFILE.write("\n")
    for fqgz_file in fqgz_list:
        print(fqgz_file)
        #basename, sample_name, silli, ind_num = parse_name_fqgz(fqgz_file)
        print(int(gzipFileSize(fqgz_file)))
        OUTFILE.write("\t".join(parse_name_fqgz(fqgz_file)))
        OUTFILE.write("\t" + str(int(gzipFileSize(fqgz_file))))
        OUTFILE.write("\t" + str(nlines_gzip(fqgz_file)/4.0))
        OUTFILE.write("\t" + str(lenreads_gzip(fqgz_file)))
        OUTFILE.write("\n")
        OUTFILE.flush()