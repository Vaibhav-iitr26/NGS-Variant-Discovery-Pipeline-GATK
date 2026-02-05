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

## Setup
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

