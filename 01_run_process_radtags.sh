process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CM01' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CM01" -b "/media/Shared/Data/chum/populations/barcodes/CM01_barcodes.txt" -y zfastq 	-e sbfI -i gzfastq -c -q -r -t 94

# notice phred64 encoding
process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CM05' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CM05" -b "/media/Shared/Data/chum/populations/barcodes/CM05_barcodes.txt" -y fastq 	-e sbfI -i gzfastq -c -q -r -t 94 -E phred64

process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CM06' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CM06" -b "/media/Shared/Data/chum/populations/barcodes/CM06_barcodes.txt" -y fastq 	-e sbfI -i gzfastq -c -q -r -t 94

process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CM09' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CM09" -b "/media/Shared/Data/chum/populations/barcodes/CM09_barcodes.txt" -y fastq 	-e sbfI -i gzfastq -c -q -r -t 94

# second batch of populations
process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CMX1' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CMX1" -b "/media/Shared/Data/chum/populations/barcodes/CMX1_barcodes.txt" -y fastq 	-e sbfI -i gzfastq -c -q -r -t 94

process_radtags -p '/media/Shared/Data/chum/populations/RawSeqs/gzip/CMX3' 	-o "/media/Shared/Data/chum/populations/cleanSeqs/CMX3" -b "/media/Shared/Data/chum/populations/barcodes/CMX3_barcodes.txt" -y fastq 	-e sbfI -i gzfastq -c -q -r -t 94




