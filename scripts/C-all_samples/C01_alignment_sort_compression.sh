#!/usr/bin/env bash

cd ~/ngs_course

SAMPLES="mother father son"

for sample in $SAMPLES
do
	bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
	data/fastq/"$sample"_R1.fastq.gz \
	data/fastq/"$sample"_R2.fastq.gz \
	| samtools sort \
	| samtools view -bh > results/alignments/"$sample".bam
done
