#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -N BIOS60132_Practical_Five

#Load necessary module
module load bio/2.0

#Set input path
inPath="../Data/backup"

#Run fastqc on all untrimmed files
fastqc -t 8 "$inPath"/untrimmed_fastq/*.gz

#Make trimming output directory
mkdir trimmed_fastq

#Loop over untrimmed files
echo "Beginning trimming..."
for f in "$inPath"/untrimmed_fastq/*1.fastq.gz; do
	#Clean up file name
	trimName=$(echo "$f" | sed 's/1\.fastq\.gz//g')
	fileName=$(basename "$f" 1.fastq.gz)
	#Run trimmomatic to generate trimmed files
	trimmomatic PE -threads 8 "$f" "$trimName"2.fastq.gz \
		trimmed_fastq/"$fileName"trimmed_1.fastq.gz trimmed_fastq/"$fileName"1.fastq.gz \
		trimmed_fastq/"$fileName"trimmed_2.fastq.gz trimmed_fastq/"$fileName"2.fastq.gz \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:"$inPath"/untrimmed_fastq/NexteraPE-PE.fa:2:40:15
done

#Create indexed reference
echo "Indexing reference..."
ref="$inPath"/ref_genome/ecoli_rel606.fasta
bwa index "$ref"

#Make output directories
mkdir results results/sam results/bam results/bcf results/vcf

#Loop over trimmed files
for f in trimmed_fastq/*trimmed_1.fastq.gz; do
	#Clean up file name
	trimName=$(echo "$f" | sed 's/1\.fastq\.gz//g')
	fileName=$(basename "$f" trimmed_1.fastq.gz)
	#Run necessary tools to prepare for variant calling
	echo "Preparing $fileName for variant calling..."
	bwa mem -t 8 "$ref" "$f" "$trimName"2.fastq.gz > results/sam/"$fileName"aligned.sam
	samtools view -@ 8 -S -b results/sam/"$fileName"aligned.sam > results/bam/"$fileName"aligned.bam
	samtools sort -@ 8 -o results/bam/"$fileName"sorted.bam results/bam/"$fileName"aligned.bam
	samtools index results/bam/"$fileName"sorted.bam
	#Begin identifying SNPs
	echo "Identifying SNPs..."
	bcftools mpileup --threads 8 -O b -o results/bcf/"$fileName"raw.bcf -f "$ref" results/bam/"$fileName"sorted.bam
	bcftools call --threads 8 --ploidy 1 -m -v -o results/vcf/"$fileName"variants.vcf results/bcf/"$fileName"raw.bcf 
	vcfutils.pl varFilter results/vcf/"$fileName"variants.vcf > results/vcf/"$fileName"finalVariants.vcf
done
