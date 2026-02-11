#!/usr/bin/env bash

cd ~/ngs_course/results/variants

snpEff -Xmx4g \
-v \
-dataDir /data/ \
GRCh38.99 \
trio.filtered.vcf \
> trio.filtered.snpeff.vcf
