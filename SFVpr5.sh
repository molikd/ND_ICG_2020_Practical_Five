#!/bin/bash
#$ -M svandext@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -N SFVpr5

module load bio/2.0

#fastqc -t 8 data/untrimmed/*.gz

#trim sequences
#trimmomatic syntax: 
#for infile in data/untrimmed/*_1.fastq.gz
#do
#	base=$(basename ${infile} _1.fastq.gz)
#	trimmomatic PE data/untrimmed/${base}_1.fastq.gz data/untrimmed/${base}_2.fastq.gz \
#	data/trimmed/${base}_1.trim.fastq.gz data/trimmed/${base}_1.un.trim.fastq.gz \
#	data/trimmed/${base}_2.trim.fastq.gz data/trimmed/${base}_2.un.trim.fastq.gz \
#	SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
#done


#index the genome. bwa is the alignment program, index is the command, argument is fasta containing reference genome
#bwa index data/ref_genome/ecoli_rel606.fasta

#align sequences to reference
#note to self...bash can't tolerate spaces before and after equal signs like R can
for x in data/trimmed/*_1.trim.fastq.gz
do
	echo "on file $x"
	base=$(basename ${x} _1.trim.fastq.gz)
	echo "basename is $base"
	forward=data/trimmed/${base}_1.trim.fastq.gz
	reverse=data/trimmed/${base}_2.trim.fastq.gz
	sam=data/results/sam/${base}.aligned.sam
	bam=data/results/bam/${base}.aligned.bam
	sorted_bam=data/results/bam/${base}.aligned.sorted.bam
	raw_bcf=data/results/bcf/${base}_raw.bcf
	variants=data/results/vcf/${base}_variants.vcf
	final_variants=data/results/vcf/${base}_final_variants.vcf

	bwa mem -t 8 data/ref_genome/ecoli_rel606.fasta $forward $reverse > $sam
#convert to bam
	samtools view -@ 8 -S -b $sam > $bam
	samtools sort -@ 8 -o $sorted_bam $bam
	samtools index $sorted_bam
#calculate read position
	bcftools mpileup --threads 8 -O b -o $raw_bcf -f data/ref_genome/ecoli_rel606.fasta $sorted_bam
#call snps, ploidy = 1 because prokaryotic chromosomes haploid. generates call format files.
	bcftools call --threads 8 --ploidy 1 -m -v -o $variants $raw_bcf 
#filter snps, generate final vcf
	vcfutils.pl varFilter $variants > $final_variants
done
