clone_filter \
	-1 /media/Shared/Data/chum/PE/gz/CMUWHAP_P_s_3_1_sequence.txt.gz \
	-2 /media/Shared/Data/chum/PE/gz/CMUWHAP_P_s_3_2_sequence.txt.gz \
-i gzfastq -o /media/Shared/Data/chum/PE/clone_filtered -D


# for just P1 sequence
process_radtags \
	-f /media/Shared/Data/chum/PE/clone_filtered/CMUWHAP_P_s_3_1_sequence.txt.fil.fq_1 \
	-b /media/Shared/Data/chum/PE/barcodes/chum_PE_barcodes.txt \
	-o /media/Shared/Data/chum/populations/cleanSeqs/CMUW \
	-i gzfastq -y fastq -c -q -r -E phred64 --filter_illumina -t 94 -e sbfI

# for PE data 
process_radtags \
	-1 /media/Shared/Data/chum/PE/clone_filtered/CMUWHAP_P_s_3_1_sequence.txt.fil.fq_1 \
	-2 /media/Shared/Data/chum/PE/clone_filtered/CMUWHAP_P_s_3_2_sequence.txt.fil.fq_2 \
	-b /media/Shared/Data/chum/PE/barcodes/chum_PE_barcodes.txt \
	-o /media/Shared/Data/chum/populations/cleanSeqs/PE_Hoodsport \
	-i gzfastq -y fastq -c -q -r -E phred64 --filter_illumina -t 94 -e sbfI
