#!/bin/bash
#$ -M sweaver4@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -N sweaver4_Practical_Five

module load bio/2.0

#use fastqc on all the zipped files
fastqc -t 8 untrimmed_fastq/*.gz


#use a for loop to go through all the gz files and trim them
for infile in untrimmed_fastq/*_1.fastq.gz
 do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE untrimmed_fastq/${base}_1.fastq.gz untrimmed_fastq/${base}_2.fastq.gz \
                trimmed_fastq/${base}_1.trim.fastq.gz trimmed_fastq/${base}_1.un.trim.fastq.gz \
                trimmed_fastq/${base}_2.trim.fastq.gz trimmed_fastq/${base}_2.un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
 done

#index reference genome
bwa index ref_genome/ecoli_rel606.fasta


#for loop for variant calling (heavily adapted from data carpentry tutorial)
for infile in trimmed_fastq/*_1.trim.fastq.gz
	do
	echo "On file $infile"
	
	base=$(basename $infile _1.trim.fastq.gz)
	echo "basename is $base"
	
	#I already made the directories "results/sam", "results/bam", etc

	forward=trimmed_fastq/${base}_1.trim.fastq.gz
	reverse=trimmed_fastq/${base}_2.trim.fastq.gz
	sam=results/sam/${base}.aligned.sam
	bam=results/bam/${base}.aligned.bam
	sorted_bam=results/bam/${base}.aligned.sorted.bam
	raw_bcf=results/bcf/${base}_raw.bcf
	variants=results/vcf/${base}_variants.vcf
	final_variants=results/vcf/${base}_final_variants.vcf

	bwa mem -t 8 ref_genome/ecoli_rel606.fasta $forward $reverse > $sam
	samtools view -@ 8 -S -b $sam > $bam
  	samtools sort -@ 8 -o $sorted_bam $bam
   	samtools index $sorted_bam
        bcftools mpileup --threads 8 -O b -o $raw_bcf -f ref_genome/ecoli_rel606.fasta $sorted_bam
   	bcftools call --threads 8 --ploidy 1 -m -v -o $variants $raw_bcf 
   	vcfutils.pl varFilter $variants > $final_variants
	
	done


