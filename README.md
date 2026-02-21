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

## Variant Calling Preparation: Indexing Files
Before variant calling, all required input files must be properly indexed. GATK relies on these index files to efficiently access genomic regions and perform statistical modeling during variant discovery. If any required index is missing or inconsistent, GATK will fail or produce incorrect results. We need index files of reference genome and known variant sites (VCF) (used for recalibration).

Known variant sites are provided as compressed VCF files and must be indexed before use in base quality score recalibration.
```
gatk IndexFeatureFile --input variants/1000g_gold_standard.indels.filtered.vcf
gatk IndexFeatureFile --input variants/GCF.38.filtered.renamed.vcf
```
Refer script `A03_create_vcf_indices.sh`

Index file for reference obtained using : 
```
samtools faidx reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
gatk CreateSequenceDictionary --REFERENCE reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
```
Refer script `A04_create_fasta_index.sh`


```
#Output
data/reference/
├── Homo_sapiens.GRCh38.dna.chromosome.20.dict
├── Homo_sapiens.GRCh38.dna.chromosome.20.fa

data/variants/
├── 1000g_gold_standard.indels.filtered.vcf
├── 1000g_gold_standard.indels.filtered.vcf.idx
├── GCF.38.filtered.renamed.vcf
├── GCF.38.filtered.renamed.vcf.idx

```
## Base Quality Score Recalibration
Base Quality Score Recalibration (BQSR) is a preprocessing step that corrects systematic errors made by the sequencing machine when estimating the accuracy of each base call.
Sequencers tend to over- or under-estimate base quality scores depending on factors such as sequencing cycle and nucleotide context. BQSR models these errors and adjusts base quality scores accordingly to improve variant calling accuracy.

BQSR helps to, reduce false-positive variant calls, improve confidence in true variants, correct machine-specific and context-dependent biases and provide more accurate input for variant calling algorithms. 

The recalibration model is built using known variant sites to distinguish real variants from sequencing errors with the help of gatk `BaseRecalibrator`
```
gatk BaseRecalibrator \
--reference <reference.fa> \
--input <alignment.bam> \
--known-sites <variants1.vcf> \
--known-sites <variants2.vcf> \
--output <output.table>
```
The recalibration model is then applied to adjust base quality scores in the BAM file.
```
gatk ApplyBQSR \
  -R reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
  -I sample.sorted.RG.dedup.bam \
  --bqsr-recal-file sample.recal_data.table \
  -O sample.sorted.RG.dedup.recal.bam
```
Recalibration performed initially on mothers sample data. Refer script `B09_perform_bqsr.sh`
And then applied for all sample, refer script `C05_perform_bqsr.sh`
```
# Output
tree results/bqsr/
results/bqsr/
├── father.recal.bai
├── father.recal.bam
├── father.recal.table
├── mother.recal.bai
├── mother.recal.bam
├── mother.recal.table
├── son.recal.bai
├── son.recal.bam
└── son.recal.table
```

## Variant Calling 

Variant calling is the process of identifying genomic positions where the sequenced sample differs from the reference genome. In this workflow, variant discovery is performed using `GATK HaplotypeCaller`, which applies local de novo assembly and probabilistic modeling to accurately detect SNPs and small indels. HaplotypeCaller is designed to distinguish true biological variants from sequencing errors, making it a standard tool for germline variant discovery.

Variants were called using HaplotypeCaller in GVCF mode, which records variant and non-variant sites.
```
gatk HaplotypeCaller \
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--input results/bqsr/mother.recal.bam \
--output results/variants/mother.HC.vcf \
--intervals chr20:10018000-10220000
```
Refer Script : B10_run_haplotype_caller.sh

```
#Output

less -S results/variants/mother.HC.vcf

##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller --output results/variants/mother.HC.vcf --intervals chr20:10018000-10220000 --input results/bqsr/mother.recal.bam --reference data/reference/>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##contig=<ID=chr20,length=64444167>
##source=HaplotypeCaller
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  mother
chr20   10019348        .       A       ACT     559.02  .       AC=2;AF=1.00;AN=2;DP=14;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=25.36;SOR=5.283        GT:AD:DP:GQ:PL  1/1:0,13:13:39:573,39,0
chr20   10019469        .       C       T       281.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.579;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=12.25;ReadPosRankSum=0.15>
chr20   10019563        .       C       T       346.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=2.315;DP=35;ExcessHet=0.0000;FS=1.313;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=9.90;ReadPosRankSum=-0.421>
chr20   10019791        .       T       G       706.06  .       AC=2;AF=1.00;AN=2;DP=29;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=24.35;SOR=1.255        GT:AD:DP:GQ:PL  1/1:0,29:29:84:720,84,0
chr20   10019950        .       T       A       997.06  .       AC=2;AF=1.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.24;SOR=0.976        GT:AD:DP:GQ:PL  1/1:0,30:30:90:1011,90,0
chr20   10020046        .       G       A       573.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=1.859;DP=39;ExcessHet=0.0000;FS=2.961;MLEAC=1;MLEAF=0.500;MQ=59.25;MQRankSum=-1.252;QD=15.10;ReadPosRankSum=0.11>
chr20   10020110        .       T       A       1266.06 .       AC=2;AF=1.00;AN=2;DP=39;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=59.18;QD=32.46;SOR=0.850        GT:AD:DP:GQ:PL  1/1:0,39:39:99:1280,117,0
chr20   10020371        .       T       G       397.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=1.442;DP=35;ExcessHet=0.0000;FS=3.064;MLEAC=1;MLEAF=0.500;MQ=42.41;MQRankSum=-3.219;QD=12.05;ReadPosRankSum=-0.5>
.
.
.
```
You can check the number of variants using 
```
grep -v '^#' results/variants/mother.HC.vcf | wc -l
411
# Total 411 variants obtained
```
After variant calling and genotyping, the resulting VCF file contains detailed information about each variant, including genomic position, reference and alternate alleles, quality metrics, and genotypes. `VariantsToTable` is used to extract selected fields from the VCF file and convert them into a tabular format for easier inspection, filtering, and downstream analysis.
```
gatk VariantsToTable \
--variant variants/mother.HC.vcf \
--fields CHROM -F POS -F TYPE -GF GT \
--output variants/mother.HC.table
```
Refer Script : B11_variants_to_table.sh

```
#Output

CHROM   POS     TYPE    mother.GT
chr20   10019348        INDEL   ACT/ACT
chr20   10019469        SNP     C/T
chr20   10019563        SNP     C/T
chr20   10019791        SNP     G/G
chr20   10019950        SNP     A/A
chr20   10020046        SNP     G/A
chr20   10020110        SNP     A/A
chr20   10020371        SNP     T/G
chr20   10020650        SNP     A/A
chr20   10020788        INDEL   AAGGCT/AAGGCT
chr20   10020826        SNP     T/T
.
.
.
.
.
.
```
```
cut -f 3 results/variants/mother.HC.table | tail -n +2 | sort | uniq -c
     84 INDEL
      1 MIXED
    326 SNP
```

This is variant calling. Now we will perform on all samples. 
Here we called variants in GVCF mode. `--emit-ref-confidence GVCF` enables GVCF mode, which records variant and non-variant sites. Later we want to combine the variant calls. For efficient merging of vcfs, we will need to output the variants as a GVCF. The `--bam-output `option in GATK HaplotypeCaller generates an additional BAM file containing reads that were locally realigned during variant calling. This file is mainly used for visualization and debugging of variant calls and is not used for downstream analysis.
```
 gatk HaplotypeCaller \
    --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    --input results/bqsr/"$SAMPLE".recal.bam \
    --output results/variants/"$SAMPLE".HC.g.vcf \
    --bam-output results/variants/"$SAMPLE".phased.bam \
    --intervals chr20:10018000-10220000 \
    --emit-ref-confidence GVCF
```
Refer Script `C06_run_haplotypecaller.sh`

## Combining GVCFs

When variant calling is performed in GVCF mode, each sample produces an intermediate GVCF file containing variant and non-variant site information. These files must be combined before final genotyping. Combining GVCFs allows multiple samples to be analyzed together using a consistent genomic framework.

Instead of directly combining GVCF files, GATK recommends importing them into a GenomicsDB workspace. GenomicsDB is a scalable storage format optimized for joint genotyping of multiple samples. This approach is more efficient and robust, especially when working with multiple samples or large genomic regions.

You can generate a GenomicsDB on our samples using `GenomicsDBImport`
```
gatk GenomicsDBImport \
--variant results/variants/mother.HC.g.vcf \
--variant results/variants/father.HC.g.vcf \
--variant results/variants/son.HC.g.vcf \
--intervals chr20:10018000-10220000 \
--genomicsdb-workspace-path results/genomicsdb
```
Refer Script `C07_create_genomicsdb.sh`

You can retrieve the combined vcf from the database with `gatk GenotypeGVCFs`
```
atk GenotypeGVCFs \
--reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--variant gendb://results/genomicsdb \
--intervals chr20:10018000-10220000 \
--output results/variants/trio.vcf
```
Refer Script `C08_genotype_gvcfs.sh`

```
#Output
 less -S results/variants/trio.vcf
.
.
.
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  father  mother  son
chr20   10019252        .       G       C       134.68  .       AC=1;AF=0.167;AN=6;BaseQRankSum=1.29;DP=15;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=16.83;ReadPosRankSum=0.439;S>
chr20   10019348        .       A       ACT     1587.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.160;DP=45;ExcessHet=0.0000;FS=3.837;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=25.36;ReadPosRankSum=0.413;>
chr20   10019469        .       C       T       1792.98 .       AC=4;AF=0.667;AN=6;BaseQRankSum=2.00;DP=89;ExcessHet=0.9691;FS=5.718;MLEAC=4;MLEAF=0.667;MQ=60.00;MQRankSum=0.00;QD=21.09;ReadPosRankSum=2.21;SO>
chr20   10019563        .       C       T       1905.98 .       AC=4;AF=0.667;AN=6;BaseQRankSum=2.32;DP=97;ExcessHet=0.9691;FS=0.809;MLEAC=4;MLEAF=0.667;MQ=60.00;MQRankSum=0.00;QD=20.06;ReadPosRankSum=-4.210e>
chr20   10019791        .       T       G       1998.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.45;DP=100;ExcessHet=0.0000;FS=0.000;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=19.99;ReadPosRankSum=-7.190>
chr20   10019950        .       T       A       3279.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=-6.280e-01;DP=108;ExcessHet=0.0000;FS=0.000;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=30.37;ReadPosRankSum=>
chr20   10020046        .       G       A       2521.98 .       AC=4;AF=0.667;AN=6;BaseQRankSum=2.21;DP=111;ExcessHet=0.9691;FS=0.000;MLEAC=4;MLEAF=0.667;MQ=57.83;MQRankSum=-1.252e+00;QD=23.35;ReadPosRankSum=>
chr20   10020110        .       T       A       3300.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.81;DP=121;ExcessHet=0.0000;FS=0.000;MLEAC=5;MLEAF=0.833;MQ=58.17;MQRankSum=-4.640e-01;QD=27.74;ReadPosRankSum=>
chr20   10020371        .       T       G       1814.98 .       AC=4;AF=0.667;AN=6;BaseQRankSum=1.44;DP=98;ExcessHet=0.9691;FS=1.741;MLEAC=4;MLEAF=0.667;MQ=41.52;MQRankSum=-3.219e+00;QD=19.31;ReadPosRankSum=->
chr20   10020650        .       T       A       2834.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=-4.260e-01;DP=106;ExcessHet=0.0000;FS=0.000;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=28.64;ReadPosRankSum=>
chr20   10020788        .       A       AAGGCT  3908.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.125;DP=113;ExcessHet=0.0000;FS=1.075;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=28.73;ReadPosRankSum=-1.68>
chr20   10020826        .       C       T       2733.93 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.41;DP=98;ExcessHet=0.0000;FS=4.011;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=28.18;ReadPosRankSum=0.397;S>
.
.
.
```
## Visualization
fter generating the final combined VCF using GenotypeGVCFs, variant calls were manually inspected using IGV (Integrative Genomics Viewer). Visualization helps verify that called variants are supported by aligned sequencing reads and are not artifacts of alignment or sequencing errors.

Following file uploaded in IGV : 

variants/mother.phased.bam

variants/mother.phased.bai

bqsr/mother.recal.bam

bqsr/mother.recal.bai

variants/mother.HC.vcf

Visualization focused in region of chr20:10,026,397-10,026,638

1.

<img width="1919" height="1031" alt="Screenshot 2026-02-21 223756" src="https://github.com/user-attachments/assets/c5ffbb59-f59c-4d30-b921-78b481f4a8ce" />

This IGV screenshot shows a region on chromosome 20 (GRCh38) after joint genotyping using GenotypeGVCFs. The VCF track highlights SNPs identified in the final variant call set. Both the recalibrated BAM and the haplotype-aware phased BAM show consistent read support for the alternate alleles at these positions. The agreement between raw alignments and haplotype-resolved reads confirms the reliability of the called variants. No obvious strand bias or alignment artifacts are observed. This visual inspection validates the accuracy of the variant calling pipeline.

2.

<img width="1919" height="1028" alt="Screenshot 2026-02-21 223922" src="https://github.com/user-attachments/assets/2f1faf20-6df1-4b29-968e-3e0d5cc8b0ea" />

This IGV snapshot shows a broader region on chromosome 20 (GRCh38) containing multiple variant sites identified after joint genotyping with GenotypeGVCFs. The VCF track displays a dense set of SNPs across the region, indicating substantial sequence variation. Both the recalibrated BAM and the phased BAM show consistent read coverage and alternate allele support at these positions. The phased BAM highlights how reads are grouped into haplotypes in variant-rich regions. Overall read depth is uniform, with no major dropouts or alignment anomalies. This visualization confirms that the variant calls are supported across a larger genomic interval.

3.

<img width="1918" height="1030" alt="Screenshot 2026-02-21 224527" src="https://github.com/user-attachments/assets/1375614e-bc07-43df-ae4c-ab5c0088d943" />

In the track with mother.phased.bam, right click on the reads and select Group alignments by > read group. This splits your track in two parts, one with artificial reads describing haplotypes that were taken in consideration (ArtificalHaplotypeRG), and one with original reads that support the variants.

Now colour the reads by phase. Do that with by right clicking on the track and choose Colour alignments by > tag, and type in “HC” (that’s the tag where the phasing is stored). When the reads are colored by the HC tag (haplotype tag), one artificial haplotype read group shows no support from the original sequencing reads, which is actually representing homozygous reference allele. The remaining reads show consistent support for the alternate alleles of both SNPs on the same haplotype, meaning the two SNPs are in phase. This confirms that the variants are correctly phased and supported by the underlying sequencing data.

