#!/bin/bash
#$ -M nvincen2@nd.edu
#$ -m abe
#$ -pe smp 4
#$ -N BIOS60132_Practical_Five

module load bio/2.0

fastqc -t 8 data/untrimmed_fastq/*.gz

for infile in data/untrimmed_fastq/*_1.fastq.gz
 do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE data/untrimmed_fastq/${base}_1.fastq.gz data/untrimmed_fastq/${base}_2.fastq.gz \
                data/trimmed_fastq/${base}_1.trim.fastq.gz data/trimmed_fastq/${base}_1.un.trim.fastq.gz \
                data/trimmed_fastq/${base}_2.trim.fastq.gz data/trimmed_fastq/${base}_2.un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
 done

bwa index data/ref_genome/ecoli_rel606.fasta

mkdir results results/sam results/bam results/bcf results/vcf

for f in ./data/trimmed_fastq/*_1.trim.fastq
do
echo "in ${f}"

base=$(basename data/trimmed_fastq/${f} _1.trim.fastq)
echo "basename is $base"
	
r1=data/trimmed_fastq/${base}_1.trim.fastq
r2=data/trimmed_fastq/${base}_2.trim.fastq
sam=results/sam/${base}.aligned.sam
bam=results/bam/${base}.aligned.bam
sorted_bam=results/bam/${base}.aligned.sorted.bam
raw_bcf=results/bcf/${base}_raw.bcf
variants=results/vcf/${base}_variants.vcf
final_variants=results/vcf/${base}_final_variants.vcf

bwa mem -t 8 data/ref_genome/ecoli_rel606.fasta $r1 $r2 > $sam
samtools view -@ 8 -S -b $sam > $bam
samtools sort -@ 8 -o $sorted_bam $bam
samtools index $sorted_bam
bcftools mpileup --threads 8 -O b -o $raw_bcf -f data/ref_genome/ecoli_rel606.fasta $sorted_bam
bcftools call --threads 8 --ploidy 1 -m -v -o $variants $raw_bcf 
vcfutils.pl varFilter $variants > $final_variants

done
