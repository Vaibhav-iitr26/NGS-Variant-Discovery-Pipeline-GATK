#!/usr/bin/env bash

cd ~/ngs_course/results

gatk SelectVariants \
--variant variants/trio.vcf \
--select-type-to-include INDEL \
--output variants/trio.INDEL.vcf

