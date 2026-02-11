#!/usr/bin/env bash

cd ~/ngs_course/results/variants

gatk SelectVariants \
--variant trio.filtered.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--output mother.trio.filtered.vcf


