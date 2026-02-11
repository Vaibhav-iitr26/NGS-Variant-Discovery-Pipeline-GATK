#!/usr/bin/env bash

cd ~/ngs_course/results/variants

gatk MergeVcfs \
--INPUT trio.SNP.filtered.vcf \
--INPUT trio.INDEL.filtered.vcf \
--OUTPUT trio.filtered.vcf

