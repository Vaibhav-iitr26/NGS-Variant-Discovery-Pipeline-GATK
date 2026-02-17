# NGS-Variant-Discovery-Pipeline-GATK
### Note - Readme inprogress

## Variant Analysis
This repository demonstrates a real-world, end-to-end NGS variant analysis pipeline, starting from raw sequencing reads to align, preprocess, discover variants, filter quality, annotate functionality, and validate results.

Variant analysis is the process of identifying genetic differences between a sequenced sample and a reference genome. These differences, called variants, mainly observed in single nucleotide changes (SNPs) and small insertions and deletions (indels).
Such variants are responsible for normal genetic diversity as well as many diseases, including cancer and inherited disorders.

In this training workflow, raw sequencing reads are processed step-by-step to produce high-confidence genetic variants.

Using the provided sequencing dataset, the pipeline:
1. Aligns raw DNA reads to a reference genome
2. Removes technical artifacts such as PCR duplicates
3. Corrects systematic sequencing errors
4. Calls SNPs and indels using probabilistic models
5. Filters low-quality variants
6. Annotates variants to understand their biological impact

The final output is a curated set of variants that represent true genetic differences in the sample.


<img width="858" height="387" alt="Screenshot 2026-02-05 175719" src="https://github.com/user-attachments/assets/a8c75966-d128-4b3b-a026-1e2a27db069d" />

---
This specific workflow and related material provided by the SIB Swiss Institute of Bioinformatics ( [NGS Variant Analysis](https://sib-swiss.github.io/NGS-variants-training/2024.9/) ), and all steps are executed in a reproducible manner using scripted workflows and controlled software environments.

## Setup and Environment Preparation
Complete workflow performed on linux based system. All required bioinformatics tools were installed using conda. All analyses were performed in a controlled software environment to ensure reproducibility.

Created `environment.yml` file
  ```
  cat environment.yml
  name: ngs-tools
  channels:
    - defaults
  dependencies:
    - samtools=1.12
    - bwa=0.7.17
    - r-base
    - snpeff=5.0
    - gatk4=4.2.6
    - python=3.8
  ```

Generated the conda environment
```
conda env create --name ngs-tools -f environment.yml
```

Activated the environment
```
conda activate ngs-tools
```
All the required directories created using `mkdir`
```
mkdir script
cd script
  mkdir A-prepare_references B-mother_only C-all_samples
```

<!--
├── data

│   ├── fastq

│   │   ├── father_R1.fastq.gz

│   │   ├── father_R2.fastq.gz

│   │   ├── mother_R1.fastq.gz

│   │   ├── mother_R2.fastq.gz

│   │   ├── son_R1.fastq.gz

│   │   └── son_R2.fastq.gz

│   ├── reference

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.dict

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa.amb

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa.ann

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa.bwt

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa.fai

│   │   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa.pac

│   │   └── Homo_sapiens.GRCh38.dna.chromosome.20.fa.sa

│   └── variants

│       ├── 1000g_gold_standard.indels.filtered.vcf

│       ├── 1000g_gold_standard.indels.filtered.vcf.idx

│       ├── GCF.38.filtered.renamed.vcf

│       ├── GCF.38.filtered.renamed.vcf.idx

│       ├── NA12878.vcf.gz

│       └── NA12878.vcf.gz.tbi

├── directories.txt

├── environment.yml

├── results

│   ├── alignments

│   │   ├── father.bam

│   │   ├── father.rg.bam

│   │   ├── father.rg.md.bam

│   │   ├── father.rg.md.bam.bai

│   │   ├── marked_dup_metrics_father.txt

│   │   ├── marked_dup_metrics_mother.txt

│   │   ├── marked_dup_metrics_son.txt

│   │   ├── mother.bam

│   │   ├── mother.rg.bam

│   │   ├── mother.rg.md.bam.flagstat

│   │   ├── mother.sam

│   │   ├── mother.sam.flagstat

│   │   ├── mother.sorted.sam

│   │   ├── son.bam

│   │   ├── son.rg.bam

│   │   ├── son.rg.md.bam

│   │   └── son.rg.md.bam.bai

│   ├── bqsr

│   │   ├── father.recal.bai

│   │   ├── father.recal.bam

│   │   ├── father.recal.table

│   │   ├── mother.recal.bai

│   │   ├── mother.recal.bam

│   │   ├── mother.recal.table

│   │   ├── son.recal.bai
│   │   ├── son.recal.bam
│   │   └── son.recal.table
│   ├── genomicsdb
│   │   ├── __tiledb_workspace.tdb
│   │   ├── callset.json
│   │   ├── chr20$10018000$10220000
│   │   │   ├── __668aa7cb-6f40-48d8-95b8-1a83938f2707133691370100416_1769869305299
│   │   │   │   ├── AD.tdb
│   │   │   │   ├── AD_var.tdb
│   │   │   │   ├── ALT.tdb
│   │   │   │   ├── ALT_var.tdb
│   │   │   │   ├── BaseQRankSum.tdb
│   │   │   │   ├── DP.tdb
│   │   │   │   ├── DP_FORMAT.tdb
│   │   │   │   ├── END.tdb
│   │   │   │   ├── ExcessHet.tdb
│   │   │   │   ├── FILTER.tdb
│   │   │   │   ├── FILTER_var.tdb
│   │   │   │   ├── GQ.tdb
│   │   │   │   ├── GT.tdb
│   │   │   │   ├── GT_var.tdb
│   │   │   │   ├── ID.tdb
│   │   │   │   ├── ID_var.tdb
│   │   │   │   ├── InbreedingCoeff.tdb
│   │   │   │   ├── MIN_DP.tdb
│   │   │   │   ├── MLEAC.tdb
│   │   │   │   ├── MLEAC_var.tdb
│   │   │   │   ├── MLEAF.tdb
│   │   │   │   ├── MLEAF_var.tdb
│   │   │   │   ├── MQRankSum.tdb
│   │   │   │   ├── PGT.tdb
│   │   │   │   ├── PGT_var.tdb
│   │   │   │   ├── PID.tdb
│   │   │   │   ├── PID_var.tdb
│   │   │   │   ├── PL.tdb
│   │   │   │   ├── PL_var.tdb
│   │   │   │   ├── PS.tdb
│   │   │   │   ├── QUAL.tdb
│   │   │   │   ├── RAW_MQandDP.tdb
│   │   │   │   ├── REF.tdb
│   │   │   │   ├── REF_var.tdb
│   │   │   │   ├── ReadPosRankSum.tdb
│   │   │   │   ├── SB.tdb
│   │   │   │   ├── __book_keeping.tdb.gz
│   │   │   │   ├── __coords.tdb
│   │   │   │   └── __tiledb_fragment.tdb
│   │   │   ├── __array_schema.tdb
│   │   │   └── genomicsdb_meta_dir
│   │   │       ├── genomicsdb_column_bounds.json
│   │   │       └── genomicsdb_meta_143e26b4-9c53-4b43-9d1d-36c57e9fdf4a.json
│   │   ├── vcfheader.vcf
│   │   └── vidmap.json
│   ├── sample_rg_fields.txt
│   └── variants
│       ├── concordance.mother.trio
│       ├── concordance.mother.trio.filtered
│       ├── father.HC.g.vcf
│       ├── father.HC.g.vcf.idx
│       ├── father.phased.bai
│       ├── father.phased.bam
│       ├── mother.HC.g.vcf
│       ├── mother.HC.g.vcf.idx
│       ├── mother.HC.table
│       ├── mother.HC.vcf
│       ├── mother.HC.vcf.idx
│       ├── mother.phased.bai
│       ├── mother.phased.bam
│       ├── mother.trio.filtered.vcf
│       ├── mother.trio.filtered.vcf.idx
│       ├── mother.trio.vcf
│       ├── mother.trio.vcf.idx
│       ├── snpEff_genes.txt
│       ├── snpEff_summary.html
│       ├── son.HC.g.vcf
│       ├── son.HC.g.vcf.idx
│       ├── son.phased.bai
│       ├── son.phased.bam
│       ├── trio.INDEL.filtered.vcf
│       ├── trio.INDEL.filtered.vcf.idx
│       ├── trio.INDEL.vcf
│       ├── trio.INDEL.vcf.idx
│       ├── trio.SNP.filtered.vcf
│       ├── trio.SNP.filtered.vcf.idx
│       ├── trio.SNP.vcf
│       ├── trio.SNP.vcf.idx
│       ├── trio.filtered.snpeff.vcf
│       ├── trio.filtered.vcf
│       ├── trio.filtered.vcf.idx
│       ├── trio.vcf
│       └── trio.vcf.idx
└── scripts
    ├── A-prepare_references
    │   ├── A01_download_course_data.sh
    │   ├── A02_create_bwa_index.sh
    │   ├── A03_create_vcf_indices.sh
    │   └── A04_create_fasta_index.sh
    ├── B-mother_only
    │   ├── B01_alignment.sh
    │   ├── B02_get_alignment_statistics.sh
    │   ├── B03_sort_alignment.sh
    │   ├── B04_compress_alignment.sh
    │   ├── B05_add_readgroups.sh
    │   ├── B06_mark_duplicates.sh
    │   ├── B07_get_alignment_stats_after_md.sh
    │   ├── B08_index_alignment.sh
    │   ├── B09_perform_bqsr.sh
    │   ├── B10_run_haplotype_caller.sh
    │   └── B11_variants_to_table.sh
    ├── C-all_samples
    │   ├── C01_alignment_sort_compression.sh
    │   ├── C02_add_readgroups.sh
    │   ├── C03_mark_duplicates_index.sh
    │   ├── C05_perform_bqsr.sh
    │   ├── C06_run_haplotypecaller.sh
    │   ├── C07_create_genomicdb.sh
    │   ├── C08_genotype_gvcfs.sh
    │   ├── C09_select_SNPs.sh
    │   ├── C10_select_INDELs.sh
    │   ├── C11_filter_SNPs.sh
    │   ├── C12_filter_INDELs.sh
    │   ├── C13_merge_filtered.sh
    │   ├── C14_extract_mother_only.sh
    │   ├── C15_evaluate_concordance.sh
    │   ├── C16_extract_mother_before_filtering.sh
    │   ├── C17_evaluate_concordance_before_filtering.sh
    │   └── C18_annotate_snpEff.sh
    └── calculate_genotype_likelihoods.R
-->

## Input Data
In this analysis, the input data used, provided in training itself, as part of the SIB NGS Variant Analysis training workflow. The dataset consists of raw paired-end sequencing reads (FASTQ format) along with the corresponding reference genome required for alignment and variant discovery.

Data downloaded using wget command
```
wget https://ngs-variants-training.s3.eu-central-1.amazonaws.com/ngs-variants-training.tar.gz
tar -xvf ngs-variants-training.tar.gz
```
Required commands stored in script `A01_download_course_data.sh`

Data included: 
1. Raw sequencing reads (.fastq / .fastq.gz)
2. Reference genome (.fa)
3. Known variant resources used for preprocessing (.vcf)

Due to the file size limitation, these datasets not uploaded on repository. You can obtain full datasets from the original training source: https://sib-swiss.github.io/NGS-variants-training

## Indexing Reference genome
Before read alignment and variant calling, the reference genome must be indexed. Indexing allows bioinformatics tools to quickly access specific genomic regions without scanning the entire genome file. And that's make it time and memory efficient.

Software bwa, used for indexing.
```
bwa index Homo_sapiens.GRCh38.dna.chromosome.20.fa
```
Commands stored in scripts `A02_create_bwa_index.sh`

## Read Alignment
Read alignment is the process of mapping short sequencing reads to a reference genome to determine their genomic origin. Since next-generation sequencing produces millions of short DNA fragments, accurate alignment is essential for reliable variant detection. In this workflow, paired-end sequencing reads are aligned to the reference genome using BWA-MEM, a widely used algorithm optimized for high-throughput short reads. So it will generate an alignment file that can be used for downstream variant analysis.

Paired-end FASTQ files were aligned to the indexed reference genome as follows:
```
bwa mem \
  reference.fa \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz > sample.sam
```
First we applied this on mother samples only. The respective script stored in `B01_alignment.sh` and the results stored in folder `results/alignments/`.

## Alignment Statistics
After read alignment, basic alignment statistics were generated to assess the quality and reliability of the mapping. These metrics help verify that the sequencing reads were correctly aligned to the reference genome before proceeding to variant calling. Evaluating alignment statistics at this stage is critical, as poor alignment quality directly affects downstream variant detection.

Alignment statistics were generated using `samtools flagstat`
```
samtools flagstat mother.sam > mother.sam.flagstat
```
Script > `B02_get_alignment_statistics.sh`, Results in `results/alignments/` 


