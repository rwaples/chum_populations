#GATK 

# use picard to CreateSequenceDictionary
java -jar /home/ipseg/Programs/picard/picard/dist/picard.jar CreateSequenceDictionary REFERENCE=/media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt OUTPUT=/media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt.dict

# look into how to create intervals files, and if they can be reused from freebayes
# see: http://gatkforums.broadinstitute.org/discussion/1204/what-input-files-does-the-gatk-accept-require
# likely make one interval file for ploidy == 2, and one for ploidy == 4
# could just be specified at the chr level
# -XL will exclude the intervals (i.e. blacklist)

# Look into "read-back phasing"
    # pass in the vcf + the BAM file

# using HaplotypeCaller
# here on a single sample, 

# note "fix_misencoded_quality_scores"
java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
--fix_misencoded_quality_scores \
-T HaplotypeCaller \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0008.bam \
--emitRefConfidence GVCF \
-variant_index_type LINEAR \
--variant_index_parameter 128000 \
-L /media/Shared/Data/chum/populations/fb/batch_42_CURATED_nosample.bed \
-o /home/ipseg/Desktop/waples/chum_populations/results/batch_42_GATK.raw.snps.indels.g.vcf


# three samples with UnifiedGenotyper
java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
--fix_misencoded_quality_scores \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0008.bam \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0009.bam \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0001.bam \
-o /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/GATK/UW-08-09-01_UG.vcf \
-nct 7 --read_buffer_size 1000000

# ReadBackedPhasing
# does this assume diploidy?
java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
--fix_misencoded_quality_scores \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0008.bam \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0009.bam \
-I /media/Shared/Data/chum/populations/aln/curated/CMUW10_0001.bam \
--variant /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/GATK/UW-08-09-01_UG.vcf \
--phaseQualityThresh 20 \
-maxDistMNP 100 \
--enableMergePhasedSegregatingPolymorphismsToMNP \
-o /home/ipseg/Desktop/waples/chum_populations/results/batch_42_CURATED/GATK/phased_UW-08-09-01_UG.vcf





# trying to speed up with multithreading
java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T PrintReads \
--fix_misencoded_quality_scores \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /home/ipseg/Data/chum/test/CMUW10_0008.bam \
-I /home/ipseg/Data/chum/test/CMUW10_0009.bam \
-I /home/ipseg/Data/chum/test/CMUW10_0001.bam \
-o /home/ipseg/Data/chum/test/merged.bam \
--read_filter MappingQualityZero


java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /home/ipseg/Data/chum/test/merged.bam \
-o /home/ipseg/Data/chum/test/merged.vcf \

java -jar /home/ipseg/Programs/GATK/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
-R /media/Shared/Data/chum/populations/aln/curated/ref/batch_42_CURATED.fasta.txt \
-I /home/ipseg/Data/chum/test/merged.bam \
--variant /home/ipseg/Data/chum/test/merged.vcf \
--phaseQualityThresh 20 \
-maxDistMNP 100 \
-o /home/ipseg/Data/chum/test/phased_merged.vcf

--enableMergePhasedSegregatingPolymorphismsToMNP \
