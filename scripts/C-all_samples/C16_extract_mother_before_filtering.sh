#!/usr/bin/env bash

cd ~/ngs_course/results/variants

gatk SelectVariants \
--variant trio.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--select-type-to-include SNP \
--select-type-to-include INDEL \
--output mother.trio.vcf
