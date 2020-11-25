#!/bin/bash
#$ -M msievert@nd.edu
#$ -m abe
#$ -pe smp 4
#$ -N pac5

module load bio/2.0

fastqc -t 8 data/untrimmed_fastq/*.gz
#trim fastq.gz files
for infile in data/untrimmed_fastq/*_1.fastq.gz
	do
    	base=$(basename ${infile} _1.fastq.gz)
    	trimmomatic PE untrimmed_fastq/${base}_1.fastq.gz untrimmed_fastq/${base}_2.fastq.gz \
		trimmed_fastq/${base}_1.trim.fastq.gz trimmed_fastq/${base}_1.un.trim.fastq.gz \
		trimmed_fastq/${base}_2.trim.fastq.gz trimmed_fastq/${base}_2.un.trim.fastq.gz \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
	done

#make directories for results
mkdir -p results/sam results/bam results/bcf results/vcf

#index reference genome
bwa index data/ref_genome/ecoli_rel606.fasta

for infile in data/timmed_fastq/*_1.trim.fastq.gz
    do
    echo "working with file ${infile}"

    base=$(basename ${infile} _1.trim.fastq.gz)
    echo "base name is $base"

    fq1 = data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2 = data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam = results/sam/${base}.aligned.sam
    bam = results/bam/${base}.aligned.bam
    sorted_bam = results/bam/${base}.aligned.sorted.bam
    raw_bcf = results/bcf/${base}_raw.bcf
    variants = results/bcf/${base}_variants.vcf
    final_variants=results/vcf/${base}_final_variants.vcf 

    bwa mem -t 8 data/ref_genome/ecoli_rel606.fasta $fq1 $fq2 > $sam
    samtools view -@ 8 -S -b $sam > $bam
    samtools sort -@ 8 -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup --threads 8 -O b -o $raw_bcf -f $raw_bcf -f ref_genome/ecoli_rel606.fasta $sorted_bam
    bcftools call --threads 8 --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
