#!/usr/bin/env bash

cd ~/ngs_course/data

samtools faidx reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa

gatk CreateSequenceDictionary --REFERENCE reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
