#!/bin/bash
  
if [ $# -ne 1 ];then
        echo "#usage: sh $0 [samplename]"
        exit
fi

SAMPLE=$1
BWA="/data/etc/bwa/bwa"
SAMTOOLS="/data/etc/samtools/bin/samtools"
REFERENCE="/data/reference/ucsc.hg19.fasta"
JAVA="/usr/bin/java"
PICARD="/data/etc/picard/picard.jar"
GATK="/data/etc/gatk/GenomeAnalysisTK.jar"
SNPEFF="/data/etc/snpEff/snpEff.jar"
MILLS="/data/etc/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
A1KG="/data/etc/bundle/1000G_phase1.indels.hg19.sites.vcf"
DBSNP138="/data/etc/bundle/dbsnp_138.hg19.vcf"


#=====4.Map to Reference=====#
#4.1#========================= BWA mem : FASTQ to SAM 
#고치기 bwa mem -R "@RG\tID:test\tSM:${SAMPLE}\tPL:ILLUMINA" ${REFERENCE} ${SAMPLE}_1.filt.fastq.gz ${SAMPLE}_2.filt.fastq.gz > ${SAMPLE}.mapped.sam

#4.2#========================= Samtools : SAM to BAM

#4.3#========================= Samtools sort : Make Sorted BAM

#Run from 4.1 to 4.3 in one command
#${BWA} mem -R "@RG\tID:test\tSM:${SAMPLE}\tPL:ILLUMINA" ${REFERENCE} ${SAMPLE}_1.filt.fastq.gz ${SAMPLE}_2.filt.fastq.gz | ${SAMTOOLS} view -Sb - | ${SAMTOOLS} sort - > ${SAMPLE}.sorted.bam

#=====5. Mark Duplicate======#
#5.1#=========================Picard MarkDuplicate : Sorted BAM to Markdup BAM
#${JAVA} -jar ${PICARD} MarkDuplicates I=${SAMPLE}.sorted.bam O=${SAMPLE}.markdup.bam M=${SAMPLE}.markdup.metrics.txt
#5.2#=========================Samtools index : Make BAM index
#${SAMTOOLS} index ${SAMPLE}.markdup.bam

#=====6.GATK======#
#6.1===============GATK Target Realign
#${JAVA} -jar ${GATK} -T RealignerTargetCreator -R ${REFERENCE} -I ${SAMPLE}.markdup.bam -known ${MILLS} -known ${A1KG} -o ${SAMPLE}.intervals #-L bed
#${JAVA} -jar ${GATK} -T IndelRealigner -R ${REFERENCE} -I ${SAMPLE}.markdup.bam -known ${MILLS} -known ${A1KG} -targetIntervals ${SAMPLE}.intervals -o ${SAMPLE}.realign.bam #-L bed

#6.2===============GATK BaseRecalibrator
#${JAVA} -jar ${GATK} -T BaseRecalibrator -R ${REFERENCE} -I ${SAMPLE}.realign.bam -knownSites ${MILLS} -knownSites ${A1KG} -knownSites ${DBSNP138} -o ${SAMPLE}.table #-L bed
#${JAVA} -jar ${GATK} -T PrintReads -R ${REFERENCE} -I ${SAMPLE}.realign.bam -o ${SAMPLE}.recal.bam -BQSR ${SAMPLE}.table #-L bed

#6.2===============GATK HaplotypeCaller, GenotypeGVCFs
${JAVA} -jar ${GATK} -T HaplotypeCaller -R ${REFERENCE} -I ${SAMPLE}.recal.bam --emitRefConfidence GVCF --dbsnp ${DBSNP138} -o ${SAMPLE}.g.vcf #-L bed
${JAVA} -jar ${GATK} -T GenotypeGVCFs -R ${REFERENCE} -V ${SAMPLE}.g.vcf -o ${SAMPLE}.raw.vcf #-L bed

#6.3===============GATK Variant Filter
#### Select SNP
#${JAVA} -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${SAMPLE}.raw.vcf -o ${SAMPLE}.raw.snp.vcf --selectTypeToInclude SNP

#### Filter SNP
#${JAVA} -jar ${GATK} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLE}.raw.snp.vcf -o ${SAMPLE}.filtered.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MQRankSum < -12.5  ||  ReadPosRankSum < -8.0" --filterName SNP_FILTER

#### Select INDEL
#${JAVA} -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${SAMPLE}.raw.vcf -o ${SAMPLE}.raw.indel.vcf --selectTypeToInclude INDEL

#### Filter INDEL
#${JAVA} -jar ${GATK} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLE}.raw.indel.vcf -o ${SAMPLE}.filtered.indel.vcf --filterExpression "QD < 2.0  || FS > 200.0 || ReadPosRankSum < -20.0" --filterName INDEL_FILTER

#### Combine SNPs and INDELs
#${JAVA} -jar ${GATK} -T CombineVariants -R ${REFERENCE} -V ${SAMPLE}.filtered.snp.vcf -V ${SAMPLE}.filtered.indel.vcf -o ${SAMPLE}.filtered.vcf -genotypeMergeOptions UNSORTED

#### Detailed Filter options are here. https://software.broadinstitute.org/gatk/documentation/article?id=2806













