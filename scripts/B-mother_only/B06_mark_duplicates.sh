#!/usr/bin/env bash


cd ~/ngs_course/results


gatk MarkDuplicates \
--INPUT alignments/mother.rg.bam \
--OUTPUT alignments/mother.rg.md.bam \
--METRICS_FILE alignments/mark_dup_metrics_mother.txt
