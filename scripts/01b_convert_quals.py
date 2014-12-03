import glob
import os.path
import HTSeq

# converting PE data to phred33
files_to_covert = glob.glob(os.path.join("/media/Shared/Data/chum/populations/cleanSeqs/", "CMUW10_*.fq.gz"))

for xfile in files_to_covert:
    xseqs = HTSeq.FastqReader(xfile, qual_scale = 'solexa')
    with open(os.path.join("/media/Shared/Data/chum/populations/cleanSeqs/CMUW_phred33", os.path.basename(xfile)), 'w') as OUTFILE:
        for seq in xseqs:
            seq.write_to_fastq_file(OUTFILE)
          
              
files_to_covert = glob.glob( "/media/Shared/Data/chum/populations/cleanSeqs/PE_Hoodsport/CMUW10_0001.1.fq.gz")

for xfile in files_to_covert:
    xseqs = HTSeq.FastqReader(xfile, qual_scale = 'solexa')
    with open(os.path.join("/media/Shared/Data/chum/populations/cleanSeqs/CMUW_phred33", os.path.basename(xfile)), 'w') as OUTFILE:
        for seq in xseqs:
            seq.write_to_fastq_file(OUTFILE)

