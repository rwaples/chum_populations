clone_filter \
	-1 /media/Shared/Data/chum/PE/gz/CMUWHAP_P_s_3_1_sequence.txt.gz \
	-2 /media/Shared/Data/chum/PE/gz/CMUWHAP_P_s_3_2_sequence.txt.gz \
-i gzfastq -o /media/Shared/Data/chum/PE/clone_filtered -D


process_radtags \
	-1 /media/Shared/Data/chum/PE/clone_filtered/CMUWHAP_P_s_3_1_sequence.txt.fil.fq_1 \
	-2 /media/Shared/Data/chum/PE/clone_filtered/CMUWHAP_P_s_3_2_sequence.txt.fil.fq_2 \
	-b /media/Shared/Data/chum/PE/barcodes/chum_PE_barcodes.txt \
	-o /media/Shared/Data/chum/PE/clean_demultiplex \
	-i gzfastq -y gzfastq -c -q -r -E phred33 --filter_illumina --disable rad_check -e sbfI
