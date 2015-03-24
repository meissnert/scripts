#!/bin/bash

module load star
module load gatk
module load samtools
module load picard

sample=817-1401A-b_S1

fq1=/data/s3/averafastq/Geparquinto/817-1401A-b_S1_1.fastq.gz
fq2=/data/s3/averafastq/Geparquinto/817-1401A-b_S1_2.fastq.gz

RGID=12345
RGLB=RNAACCESS
RGPL=ILLUMINA
RGSM=817-1401A-b_S1

refgtf=/data/database/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
refFasta=/data/database/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

genomeDir=/data/database/Homo_sapiens/UCSC/hg19/star_genome
outDir=/data/storage/broad_test
threads=32

mills=/data/database/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf
g1000=/data/database/GATK/hg19/1000G_phase1.indels.hg19.vcf
dbsnp=/data/database/GATK/hg19/dbsnp_137.hg19.vcf

TMPDIR=/mnt/tobias

# 1ST PASS ALIGNMENT
mkdir -p $outDir/$sample/1pass

STAR \
	--genomeDir $genomeDir \
	--readFilesIn $fq1 $fq2 \
	--readFilesCommand zcat \
	--sjdbGTFfile $refgtf \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattrRGline ID:$RGID LB:$RGLB PL:$RGPL SM:$RGSM \
	--outReadsUnmapped Fastx \
	--runThreadN $threads \
	--outFileNamePrefix $outDir/$sample/1pass/

# 2ND PASS ALIGNMENT
secondPassDir=$outDir/$sample/hg19_2pass
mkdir -p $secondPassDir

STAR \
	--runMode genomeGenerate \
	--genomeDir $secondPassDir \
	--genomeFastaFiles $refFasta \
    --sjdbFileChrStartEnd $outDir/$sample/1pass/SJ.out.tab \
    --sjdbOverhang 100 \
    --runThreadN $threads

mkdir -p $outDir/$sample/2pass

STAR \
	--genomeDir $secondPassDir \
	--readFilesIn $fq1 $fq2 \
	--readFilesCommand zcat \
	--sjdbGTFfile $refgtf \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattrRGline ID:$RGID LB:$RGLB PL:$RGPL SM:$RGSM \
	--outReadsUnmapped Fastx \
	--runThreadN $threads \
	--outFileNamePrefix $outDir/$sample/2pass/

# MARK DUPLICATES
java -jar `which MarkDuplicates.jar` \
	I=$outDir/$sample/2pass/Aligned.sortedByCoord.out.bam \
	O=$outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.bam  \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=$outDir/$sample/2pass/output.metrics

# SPLIT'N'TRIM + MAPPING Q
java -jar `which GenomeAnalysisTK.jar` \
	-T SplitNCigarReads \
	-R $refFasta \
	-I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.bam \
	-o $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS \
	-fixNDN

# INDEL REALIGNMENT
java -jar `which GenomeAnalysisTK.jar` \
	-T RealignerTargetCreator \
	-R $refFasta \
	-I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.bam \
	-o $outDir/$sample/2pass/output.intervals \
	-known $mills \
	-known $g1000 \
	-nt 8

java -Djava.io.tmpdir=$TMPDIR -jar `which GenomeAnalysisTK.jar` \
	-I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.bam \
	-R $refFasta \
	-T IndelRealigner \
	-targetIntervals $outDir/$sample/2pass/output.intervals \
	-o $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
	-known $mills \
	-known $g1000 \
	--consensusDeterminationModel KNOWNS_ONLY \
	-LOD 0.4

# BQSR
java -Xmx16g -jar `which GenomeAnalysisTK.jar` \
    -T BaseRecalibrator \
    -R $refFasta \
    -I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -o $outDir/$sample/2pass/recal_data.table \
    -nct $threads

java -jar `which GenomeAnalysisTK.jar` \
    -T PrintReads \
    -R $refFasta \
    -I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -BQSR $outDir/$sample/2pass/recal_data.table \
    -o $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam

# VARIANT CALLING
java -jar `which GenomeAnalysisTK.jar` \
	-T HaplotypeCaller \
	-R $refFasta \
	-I $outDir/$sample/2pass/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-stand_emit_conf 20.0 \
	-o $outDir/$sample/2pass/${sample}_raw.vcf

# VARIANT FILTERING
java -jar `which GenomeAnalysisTK.jar` \
	-T VariantFiltration \
	-R $refFasta \
	-V $outDir/$sample/2pass/${sample}_raw.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS \
	-filter "FS > 30.0" \
	-filterName QD \
	-filter "QD < 2.0" \
	-o $outDir/$sample/2pass/${sample}_filt.vcf

exit 0
