#!/usr/bin/env bash

cd ~/ngs_course/results

samtools view -bh alignments/mother.sorted.sam > alignments/mother.bam
