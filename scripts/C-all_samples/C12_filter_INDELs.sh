#!/usr/bin/env bash

cd ~/ngs_course/results

gatk VariantFiltration \
--variant variants/trio.INDEL.vcf \
--filter-expression "QD < 2.0"                  --filter-name "QD2" \
--filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
--filter-expression "FS > 200.0"                --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
--output variants/trio.INDEL.filtered.vcf
