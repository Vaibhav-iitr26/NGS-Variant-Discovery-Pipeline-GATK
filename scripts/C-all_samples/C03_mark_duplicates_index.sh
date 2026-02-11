#!/usr/bin/env bash

cd ~/ngs_course/results/alignments

for sample in mother father son
do
	gatk MarkDuplicates \
	--INPUT "$sample".rg.bam \
	--OUTPUT "$sample".rg.md.bam \
	--METRICS_FILE marked_dup_metrics_"$sample".txt

	samtools index "$sample".rg.md.bam
done

