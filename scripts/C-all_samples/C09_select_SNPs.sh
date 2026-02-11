#!/usr/bin/env bash

cd ~/ngs_course/results

gatk SelectVariants \
--variant variants/trio.vcf \
--select-type-to-include SNP \
--output variants/trio.SNP.vcf
