# Reproducing HDAC1 Regulation in Dendritic Cells (ATAC-seq & RNA-seq Integration)

![Nextflow](https://img.shields.io/badge/workflow-Nextflow-23aa62)
![R](https://img.shields.io/badge/R-4.x-blue)

A Nextflow/Shell/R pipeline to process raw ATAC-seq reads, perform peak calling, and integrate with RNA-seq data to identify epigenetically regulated transcriptional programs in cDC1 and cDC2 subsets.

## Dependencies (Software & Data)

### Software Tools
- FastQC
- Trim Galore
- Bowtie2
- Samtools
- Bedtools
- MACS2
- deepTools

### R Packages
- DiffBind
- DESeq2
- ChIPseeker
- clusterProfiler
- motifmatchr
- JASPAR2020
- ggplot2

### Reference Genome
- Mouse genome: **mm10/GRCm38**

### External Data
- **ATAC-seq Raw Data**: GSE266581 (SRA accessions listed in `SRR_Acc_List.txt`).
- **RNA-seq Processed Counts**: GSE266583 (`*_Raw_counts.tsv.gz`).
- **Blacklist**: ENCODE mm10 blacklist v2.

## Pipeline Steps (How to Run)
Follow these sequential steps to reproduce Figure 6 (a–f) from *"The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity"*.

### Step 1: Data Acquisition & Preprocessing
1. Download raw FASTQ (Single-End reads):
   ```bash
   prefetch --option-file SRR_Acc_List.txt
   cat SRR_Acc_List.txt | xargs -n 1 fastq-dump --gzip
   ```
2. Quality control with FastQC:
   ```bash
   fastqc raw_data/*.fastq.gz -o qc_plots/
   ```
3. Adapter trimming using Trim Galore (Nextera adapters):
   ```bash
   trim_galore --nextera --fastqc -o trimmed_reads/ raw_data/*.fastq.gz
   ```

### Step 2: Alignment & Filtering
1. Align trimmed reads to mm10 using Bowtie2:
   ```bash
   bowtie2 --very-sensitive -x /path/to/mm10/index -U trimmed_reads/*.fq.gz -S alignments/sample.sam
   ```
2. Convert, sort, and filter BAMs (remove low-quality/unmapped reads):
   ```bash
   samtools view -bS alignments/sample.sam | samtools sort -o alignments/sample.sorted.bam
   samtools view -b -q 30 alignments/sample.sorted.bam \
     | samtools view -b -o filtered_bams/sample.filtered.bam
   ```
3. Remove mitochondrial reads and blacklist regions:
   ```bash
   samtools idxstats filtered_bams/sample.filtered.bam | cut -f 1 | grep -v 'chrM' \
     | xargs samtools view -b filtered_bams/sample.filtered.bam > filtered_bams/sample.noMT.bam
   bedtools intersect -v -abam filtered_bams/sample.noMT.bam -b ENCODE_mm10_blacklist.v2.bed \
     > filtered_bams/sample.final.bam
   ```

### Step 3: Peak Calling
Call ATAC peaks with MACS2:
```bash
macs2 callpeak -t filtered_bams/sample.final.bam -f BAM -g mm --nomodel --extsize 147 \
  --keep-dup auto -q 0.01 -n peaks/sample
```

### Step 4: Differential Accessibility (DA) Analysis
Identify differential accessible regions (DARs) in R using DiffBind (EdgeR with TMM normalization):
```bash
Rscript scripts/diffbind_DA_analysis.R
```
Outputs: `DARs_*.csv`, `Figure6_a_b_Stats.pdf`.

### Step 5: Multi-Omics Integration (RNA + ATAC)
Annotate peaks and integrate with RNA-seq to highlight concordant accessibility and transcriptional changes:
```bash
Rscript scripts/reproduce_fig6_full_pipeline.R
```
- Peak annotation via ChIPseeker to nearest genes.
- RNA-seq processing with DESeq2 using `*_Raw_counts.tsv.gz`.
- Overlap calling for gain-of-accessibility with upregulated genes (and loss with downregulated genes).

Outputs: scatter plots for Figure 6c and 6e.

### Step 6: Functional & Motif Enrichment
Perform GO and motif enrichment in R:
```bash
Rscript scripts/functional_motif_enrichment.R
```
- GO enrichment: `clusterProfiler`.
- Motif enrichment: `motifmatchr` with `JASPAR2020` motifs.

### Step 7: Visualization Tracks & QC
Generate coverage tracks and TSS enrichment plots:
```bash
bamCoverage -b filtered_bams/sample.final.bam -o tracks/sample.bigWig --normalizeUsing CPM
computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 \
  -S tracks/*.bigWig -R gencode.vM25.annotation.gtf -o qc_plots/tss_matrix.gz
plotHeatmap -m qc_plots/tss_matrix.gz -out qc_plots/TSS_enrichment_heatmap.pdf
```
Outputs: BigWig tracks (for IGV visualization) and heatmaps underpinning Figure 6d and 6f.

## File Structure
```
raw_data/
trimmed_reads/
alignments/
filtered_bams/
peaks/
tracks/
qc_plots/
results/
└── figures/
```

## Usage
Run the complete workflow with the provided shell and R scripts (adjust paths to match your environment):
```bash
bash run_pipeline.sh
Rscript scripts/reproduce_fig6_full_pipeline.R
```
Ensure all dependencies and external data are available before execution.
