#!/usr/bin/env bash

cd ~/ngs_course/data

gatk IndexFeatureFile --input variants/1000g_gold_standard.indels.filtered.vcf

gatk IndexFeatureFile --input variants/GCF.38.filtered.renamed.vcf

