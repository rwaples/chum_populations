
# gzip each of the individual fq files output by process_radtags for the population data
cd /media/Shared/Data/chum/populations/cleanSeqs
find . -name '*.fq' -exec gzip {} \;
find . -name '*.fq.gz' -exec mv -i {} /media/Shared/Data/chum/populations/cleanSeqs \;
			

# also the three parents from mapping families
cd /media/Shared/Data/chum/diploids/FASTQ
find . -name '*.fq' -exec gzip {} \;
find . -name '*.fq.gz' -exec mv -i {} /media/Shared/Data/chum/populations/cleanSeqs \;


# need to concatenate the re-sequenced individuals (CMX1)
cat /media/Shared/Data/chum/populations/cleanSeqs/CMX1/CMHAMM10_0002.fq.gz /media/Shared/Data/chum/populations/cleanSeqs/CMHAMM10_0002.fq.gz > /media/Shared/Data/chum/populations/cleanSeqs/re_CMHAMM10_0002.fq.gz
cat /media/Shared/Data/chum/populations/cleanSeqs/CMX1/CMKALA03_0010.fq.gz /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0010.fq.gz > /media/Shared/Data/chum/populations/cleanSeqs/re_CMKALA03_0010.fq.gz
cat /media/Shared/Data/chum/populations/re_cleanSeqs/CMX1/CMKALA03_0037.fq.gz /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0037.fq.gz > /media/Shared/Data/chum/populations/cleanSeqs/re_CMKALA03_0037.fq.gz
cat /media/Shared/Data/chum/populations/cleanSeqs/CMX1/CMKALA03_0044.fq.gz /media/Shared/Data/chum/populations/cleanSeqs/CMKALA03_0044.fq.gz > /media/Shared/Data/chum/populations/cleanSeqs/re_CMKALA03_0044.fq.gz

# rename and move the remaining files manually
