i#!/bin/bash
#$ -pe smp 8


bwa mem -t 8 ecoli_rel606.fasta  SRR2585866_1.trim.sub.fastq SRR2585866_2.trim.sub.fastq > SRR2584866.aligned.sam

samtools view -@ 8 -S -b SRR2584866.aligned.sam > SRR2584866.aligned.bam

samtools sort -@ 8 -o SRR2584866.aligned.sorted.bam SRR2584866.aligned.bam

bcftools mpileup -threads 8 -O b -o SRR2584866_raw.bcf -f ecoli_rel606.fasta SRR2584866.aligned.sorted.bam 

bcftools call --threads 8 --ploidy 1 -m -v -o SRR2584866_variants.vcf SRR2584866_raw.bcf 

vcfutils.pl varFilter SRR2584866_variants.vcf  > SRR2584866_final_variants.vcf
