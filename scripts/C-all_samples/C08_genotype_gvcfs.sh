#!/usr/bin/env bash

cd ~/ngs_project

gatk GenotypeGVCFs \
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--variant gendb://results/genomicsdb \
--intervals chr20:10018000-10220000 \
--output results/variants/trio.vcf
