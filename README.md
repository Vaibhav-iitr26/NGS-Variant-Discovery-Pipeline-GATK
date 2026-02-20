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

Output
```
~/ngs_course/data# tree
.
├── fastq
│   ├── father_R1.fastq.gz
│   ├── father_R2.fastq.gz
│   ├── mother_R1.fastq.gz
│   ├── mother_R2.fastq.gz
│   ├── son_R1.fastq.gz
│   └── son_R2.fastq.gz
├── reference
│   ├── Homo_sapiens.GRCh38.dna.chromosome.20.fa
└── variants
    ├── 1000g_gold_standard.indels.filtered.vcf
    ├── 1000g_gold_standard.indels.filtered.vcf.idx
    ├── GCF.38.filtered.renamed.vcf
    ├── GCF.38.filtered.renamed.vcf.idx
    ├── NA12878.vcf.gz
    └── NA12878.vcf.gz.tbi

4 directories, 20 files
```

Due to the file size limitation, these datasets not uploaded on repository. You can obtain full datasets from the original training source: https://sib-swiss.github.io/NGS-variants-training

## Indexing Reference genome
Before read alignment and variant calling, the reference genome must be indexed. Indexing allows bioinformatics tools to quickly access specific genomic regions without scanning the entire genome file. And that's make it time and memory efficient.

Software bwa, used for indexing.
```
bwa index Homo_sapiens.GRCh38.dna.chromosome.20.fa
```
Commands stored in scripts `A02_create_bwa_index.sh`


Output
```
~/ngs_course# ls -ltr data/reference/

total 174140
-rw-r--r-- 1  502 staff 65518298 Nov 20  2020 Homo_sapiens.GRCh38.dna.chromosome.20.fa
-rw-r--r-- 1 root root  64444268 Jan 19 15:13 Homo_sapiens.GRCh38.dna.chromosome.20.fa.bwt
-rw-r--r-- 1 root root  16111043 Jan 19 15:13 Homo_sapiens.GRCh38.dna.chromosome.20.fa.pac
-rw-r--r-- 1 root root        89 Jan 19 15:13 Homo_sapiens.GRCh38.dna.chromosome.20.fa.ann
-rw-r--r-- 1 root root      1341 Jan 19 15:13 Homo_sapiens.GRCh38.dna.chromosome.20.fa.amb
-rw-r--r-- 1 root root  32222136 Jan 19 15:13 Homo_sapiens.GRCh38.dna.chromosome.20.fa.sa
-rw-r--r-- 1 root root        24 Jan 30 13:25 Homo_sapiens.GRCh38.dna.chromosome.20.fa.fai
-rw-r--r-- 1 root root       153 Jan 30 13:25 Homo_sapiens.GRCh38.dna.chromosome.20.dict
```

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


Output
```
~/ngs_course# less -S results/alignments/mother.sam

@SQ     SN:chr20        LN:64444167
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa data/fastq/mother_R1.fastq.gz data/fastq/mother_R2.fastq.gz
H0164ALXX140820:2:1101:2136:40460       73      chr20   15996337        60      74M77S  =       15996337        0       TCAACATAATGCTGGCCCCATAAAATGAATTTGGAAGTCTTCCCTTCTCTTCAGTTTTTTTGAAAGAGTTTGAGCAAGCTAGCAAAAA>
H0164ALXX140820:2:1101:2136:40460       133     chr20   15996337        0       *       =       15996337        0       AAAAACAGGAAAAGGGCAAGAACAACGAGGACAACAAAGCAAAAGCCGAAGTGAACCGAGCGCCCACCCACACCCCCTATAAAAAAAA>
H0164ALXX140820:2:1101:2218:71015       83      chr20   10190868        60      56S95M  =       10190687        -276    TTTTTCATCCTTTTTCTTTTTCATATTTTTATTTTTTAACTTATTATTTATCACCTATTTTTTCTCAGACTTAGGACCTCCAATATTT>
H0164ALXX140820:2:1101:2218:71015       163     chr20   10190687        60      4S35M112S       =       10190868        276     AAAATGAAAGAGGAATAGAGGGGATTTTGGAGTAAGGGAAGAAAAATAATTAAGGTAAAAGGGGAGGGTATGGGATGAGG>
H0164ALXX140820:2:1101:2451:50885       83      chr20   15991863        60      57S94M  =       15991858        -99     GACTTCTAATTTTACTTTTTTTCTTTTTTTATCTTTTCCCTTTTTTCTCTTTTTTCTTTTTATTAATATTTTAGTATCTAATGAAATT>
H0164ALXX140820:2:1101:2451:50885       163     chr20   15991858        60      8S54M89S        =       15991863        99      GAATAAGATTGGGTATTATGAATATTTTAGTATATAATGAAATTAATAATATTGGTTAAGTAACTGGAAAATAAAAACAT>
H0164ALXX140820:2:1101:2522:12015       73      chr20   10063061        60      92M59S  =       10063061        0       AAGAATTAAATGAGCAATGAGTATTATTGATGTCTTAATATATGTCACACTAAAGCAACACTTCCCAACACTTCCAAAACGTAACAAG>
H0164ALXX140820:2:1101:2522:12015       133     chr20   10063061        0       *       =       10063061        0       AAATGAGGAGAATGGGTCAAAGATAGAGAAGAATAAAGACAGAGAGAAAAGGAGGAAGGATTGTGAAGATAGTTGCACATAGGACTAA>
H0164ALXX140820:2:1101:2583:10046       121     chr20   10081167        49      72S79M  =       10081167        0       CCCTTTGAATTTTTGCTACACCTAGTATAATTTTTTTTTTTCCATAGTTCCCTGATTGTTTCAATCAACATCTCAGACATGTGTTATT>
H0164ALXX140820:2:1101:2583:10046       181     chr20   10081167        0       *       =       10081167        0       GGTGCTGGCGGGGGGGCGACCTCTTGGTGCTACCCCCCCTTTTTTTCTTTCCCCCACTTTCCATTTTCCCTTACTTATTCTTTCATTT>
H0164ALXX140820:2:1101:2756:24726       99      chr20   15962356        60      86M65S  =       15962942        645     TGCAGAAACTGAGATTATGTAAAGTGCCTGAATTTCACTAGAGTTATGAAGCTTGTAGCTCTAAATTCTGAATCACAACCCTATGAAT>
H0164ALXX140820:2:1101:2756:24726       147     chr20   15962942        60      81S59M11S       =       15962356        -645    CGGGGCGGGGGGGGCCGGCGTGTCGGGTCGCAATCACAGTTCATTTTTTCATTATCTAATGATTTTTTGAATTGCTTTCA>
.
.
.
.
```


## Alignment Statistics
After read alignment, basic alignment statistics were generated to assess the quality and reliability of the mapping. These metrics help verify that the sequencing reads were correctly aligned to the reference genome before proceeding to variant calling. Evaluating alignment statistics at this stage is critical, as poor alignment quality directly affects downstream variant detection.

Alignment statistics were generated using `samtools flagstat`
```
samtools flagstat mother.sam > mother.sam.flagstat
```
Script > `B02_get_alignment_statistics.sh`

Output
```
~/ngs_course# cat results/alignments/mother.sam.flagstat
133477 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
317 + 0 supplementary
0 + 0 duplicates
132892 + 0 mapped (99.56% : N/A)
133160 + 0 paired in sequencing
66580 + 0 read1
66580 + 0 read2
131470 + 0 properly paired (98.73% : N/A)
131990 + 0 with itself and mate mapped
585 + 0 singletons (0.44% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
Alignment statistics confirm that most sequencing reads were successfully placed onto the reference genome and the data is suitable for reliable variant calling.

## Sorting and Compression of Alignment Files
The alignment output generated by the aligner is initially in SAM format, which is a plain-text representation of read alignments. While human-readable, SAM files are large and inefficient for downstream analysis. To enable efficient storage and processing, the alignment file is converted to BAM format, compressed, and sorted by genomic coordinate. Most variant analysis tools will not run on unsorted or uncompressed alignment files.

Command `samtools sort` used for sorting and `samtools view` to convert SAM file into BAM. Corresponding scripts stored in `B03_sort_alignment.sh` and `B04_compress_alignment.sh` respectively. 


Output
```
~/ngs_course# less -S results/alignments/mother.sorted.sam

@HD     VN:1.6  SO:coordinate
@SQ     SN:chr20        LN:64444167
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa data/fastq/mother_R1.fastq.gz data/fastq/mother_R2.fastq.gz
@PG     ID:samtools     PN:samtools     PP:bwa  VN:1.12 CL:samtools sort -o alignments/mother.sorted.sam alignments/mother.sam
H0164ALXX140820:2:2107:29866:23020      129     chr20   102594  0       4S69M78S        =       15950999        15848406        TCTTTCTCTCTCTCTCTCTCTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATATATA>
H0164ALXX140820:2:2124:9911:32514       129     chr20   485061  7       21S49M81S       =       10043402        9558342 TGGGCAACAGAGTGAGACGCAGTCTCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAGAAAAACTATAAATTT>
H0164ALXX140820:2:1116:7739:21157       2179    chr20   743122  0       71H67M13H       =       10151840        9408719 TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCACACACACACACACACACACNCA     @=?:@;@?0;@?@;<?>
H0164ALXX140820:2:1124:30810:46349      2179    chr20   743122  0       59H70M22H       =       10151786        9408665 TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCACACACACACACACACACACCCACAC  ?:7?==?9@><??9@:>
H0164ALXX140820:2:2113:23228:41269      2179    chr20   743122  6       58H67M26H       =       10151829        9408708 TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCACACACACACACACACACACACA     ?:?=09@-@?@?@???>
H0164ALXX140820:2:2111:17016:7936       2195    chr20   1563232 1       43H108M =       10027661        8464323 GCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGTGGACCCCCTGAGGCTGGGAGTTTCAGACCAGTCTGATCAACATGGAGAAACCTCATCTC>
H0164ALXX140820:2:1216:11474:70311      2179    chr20   1668230 0       95H24M4I12M1I7M2D8M     =       10065338        8397109 TCTTTCTTTCTTTCTTTCTTTTTCTTTTTTTTTCTTTTTCTTTTTCTTTTTCTTTT        @-08>-@,;5@;?6<8>
H0164ALXX140820:2:2108:17727:40390      2179    chr20   2077647 0       116H30M5H       =       15904614        13826968        TTTTTTTTTTTTTTTTTTTTTTTTCTTTTT  8:?,>8:>=;8;<@@4->9;,:?#######  NM:i:0  MD:Z:30 >
H0164ALXX140820:2:2212:3557:16375       2115    chr20   2077647 4       77H31M43H       =       16060745        13983099        TTTTTTTTTTTTTTTTTTTTTTTTCTTTTTG ==<<;<;;=9==,=,5,=<9<9+,,@+>,2+ NM:i:0  MD:Z:31 >
H0164ALXX140820:2:2211:28374:44627      2161    chr20   2541647 0       27H43M81H       =       10163597        7621981 AAAAAAAAAAAAAAAAAGAAAGAAAAAAAAGAAAAGAAAAAAG     ####################?,<>0>6,50-+70?-=,+;>
H0164ALXX140820:2:2108:14946:71753      2195    chr20   2825734 52      102H30M19H      =       10030650        7204888 TCAATAAAAATAAAAAAAAAAAAAAAAAAA  :+9,48,,13++,=87,=<<<<;<<=<==<  NM:i:0  MD:Z:30 MC:Z:151>
H0164ALXX140820:2:2123:22680:31389      2131    chr20   2850472 0       97H54M  =       10085518        7234994 AAAAATTATCTGGGCATGGTGTTGTGCACCTGTAGTTCCAGCTACTTTGGAGGC  ,07:3;?8<7-16>*,>,=>4<84=?,/:8?><@8>=98?>
H0164ALXX140820:2:1116:9069:68923       2115    chr20   2868619 3       77H35M39H       =       16060670        13192052        TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTCCC     ==<;:<<<======7<=7<==>>,69>+=>+--;+     >
H0164ALXX140820:2:2104:8531:46086       177     chr20   2953910 10      93S22M36S       =       15835385        12881570        TCTTCCTCATCTTTATTAGTAGTTTGGAATACATTTGTCTTTCTTTTTTTCTGTCTAATTTTTTATATTTTTAATAAATC>
H0164ALXX140820:2:2205:32576:17272      2169    chr20   3008730 0       98H30M23H       =       3008730 0       AAAAAAAAAAAAAAACAAAAAAAAAAAAAA  ################<<;<6<8<======  NM:i:0  MD:Z:30 AS:i:30 XS:i:29 >
```

sam file compressed to bam file. You can see the size difference.
```
-rw-r--r-- 1 root root 57026170 Jan 21 16:05 mother.sorted.sam
-rw-r--r-- 1 root root 16284304 Jan 25 14:17 mother.bam
```

## Adding read groups
Read groups provide metadata about how sequencing reads were generated. They describe the sample identity, sequencing platform, and run information, allowing downstream tools to correctly model technical variation. Read groups allow tools to distinguish reads from different samples or sequencing runs, So handle multi-sample and joint genotyping workflows. GATK requires valid read group information in the BAM file. Variant calling will not work correctly without it. Even for single-sample analyses, read groups are mandatory for GATK-based workflows.

Each read group contains the following key tags:
| Tag    | Description                          |
| ------ | ------------------------------------ |
| `RGID` | Read group ID (unique identifier)    |
| `RGLB` | Library name                         |
| `RGPL` | Sequencing platform (e.g., ILLUMINA) |
| `RGPU` | Platform unit (flowcell/lane)        |
| `RGSM` | Sample name                          |

Read groups were added using `gatk AddOrReplaceReadGroups`
```
gatk AddOrReplaceReadGroups \
--INPUT alignments/mother.bam \
--OUTPUT alignments/mother.rg.bam \
--RGLB lib1 \
--RGPU H0164.2.ALXX140820 \
--RGPL ILLUMINA \
--RGSM mother \
--RGID H0164.2.mother
```
Script > `B05_add_readgroups.sh`

## Marking Duplicates

During library preparation, DNA fragments are amplified using PCR. This process can generate multiple identical copies of the same original DNA fragment. These copies are called PCR duplicates. Duplicate reads do not represent independent biological evidence. If they are not handled properly, they can artificially amplify read depth and lead to false variant calls.
Marking duplicates allows downstream tools to avoid counting the same DNA fragment multiple times, so it prevent false confidence in variant evidenceand improve accuracy of variant calling.

Duplicates are marked, not removed, so the original data remains intact. 
```
gatk MarkDuplicates \
--INPUT alignments/mother.rg.bam \
--OUTPUT alignments/mother.rg.md.bam \
--METRICS_FILE alignments/marked_dup_metrics_mother.txt 
```
Stored in script `B06_mark_duplicates.sh`.

Alignment statistics obtained of duplicate marked file using `samtools flagstat`, refer `B07_get_alignment_stats_after_md.sh`
```
~/ngs_course# cat results/alignments/mother.rg.md.bam.flagstat
133477 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
317 + 0 supplementary
17329 + 0 duplicates
132892 + 0 mapped (99.56% : N/A)
133160 + 0 paired in sequencing
66580 + 0 read1
66580 + 0 read2
131470 + 0 properly paired (98.73% : N/A)
131990 + 0 with itself and mate mapped
585 + 0 singletons (0.44% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
Here you can see the number of duplicate reads.

## Indexing bam file
After marking duplicate reads, the BAM file must be indexed to allow efficient access to specific genomic regions. Indexing creates a companion index file that enables downstream tools to quickly retrieve alignments without scanning the entire BAM file. Most variant analysis and visualization tools require an indexed BAM file as input.

The BAM file was indexed using `samtools`. Refer script > `B08_index_alignment.sh` 
```
samtools index bam-file
```
`.bai ` an index file generated.
```
mother.rg.md.bam
mother.rg.md.bam.bai
```
From fasta of a sample(we used mother sample first), we have an alignment file, which sorted, having read group and duplicate marking. Also obtained their index file. A perfect for variant call. So we will apply all these actions on other samples too, and obtain required files.  

To obtain bam file from fasta, `C01_alignment_sorting_compression.sh` script executed
```
for SAMPLE in mother father son
do
    bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    data/fastq/"$SAMPLE"_R1.fastq.gz \
    data/fastq/"$SAMPLE"_R2.fastq.gz \
    | samtools sort \
    | samtools view -bh > results/alignments/"$SAMPLE".bam
done
```
To add read groups in each bam file, automated workflow used, documented in `C02_add_readgroups.sh`. `sample_rg_fields.txt` used as a input file.
```
~/ngs_course# cat results/sample_rg_fields.txt
mother  lib1    H0164.2.ALXX140820      H0164.2.mother
father  lib2    H0164.3.ALXX140820      H0164.3.father
son     lib3    H0164.6.ALXX140820      H0164.6.son
```
```
cat sample_rg_fields.txt | while read SAMPLE LB PU ID
do
    gatk AddOrReplaceReadGroups \
    --INPUT alignments/"$SAMPLE".bam \
    --OUTPUT alignments/"$SAMPLE".rg.bam \
    --RGLB "$LB" \
    --RGPU "$PU" \
    --RGPL ILLUMINA \
    --RGSM "$SAMPLE" \
    --RGID "$ID"
done 
```

Duplicates marked by executing `C03_mark_duplicates.sh` and index file of each bam file obtained by executing `C04_index_alignments.sh`.

At the end, yo will see following files 
```
~/ngs_course# tree results/alignments/
results/alignments/
├── father.bam
├── father.rg.bam
├── father.rg.md.bam
├── father.rg.md.bam.bai
├── mark_dup_metrics_mother.txt
├── marked_dup_metrics_father.txt
├── marked_dup_metrics_mother.txt
├── marked_dup_metrics_son.txt
├── mother.bam
├── mother.rg.bam
├── mother.rg.md.bam
├── mother.rg.md.bam.bai
├── mother.rg.md.bam.flagstat
├── mother.sam
├── mother.sam.flagstat
├── mother.sorted.sam
├── son.bam
├── son.rg.bam
├── son.rg.md.bam
└── son.rg.md.bam.bai

```

