# BF528 Final Project 1: ATAC-seq Analysis Pipeline

This repository contains a reproducible **Nextflow DSL2** workflow for ATAC-seq analysis of WT and HDAC1-KO cDC2 samples.

The pipeline is designed for **Project 1** requirements and includes:
- Raw read QC and report aggregation
- Optional adapter trimming
- mm10 alignment and BAM cleanup (mtDNA + blacklist + duplicates)
- ATAC-seq peak calling
- Advanced ATAC-specific QC metrics

## Project layout

- `main.nf`: Pipeline entry point
- `nextflow.config`: Parameters, profiles, and resources
- `modules/`: Tool-level process modules
- `notebooks/`: Downstream analysis report and helper code
- `environment.yml`: Conda environment for reproducibility
- `data/`, `ref/`, `results/`: Input, reference, and outputs

## Requirements

- Nextflow (>= 23)
- Conda / Mamba
- Optional: Singularity
- Linux or macOS

## Data and reference preparation

1. Download ATAC-seq FASTQ from GEO accession **GSE266584** into `data/raw_fastq/`.
2. Prepare mm10 reference files:
   - Bowtie2 index prefix: `ref/bowtie2_index/mm10`
   - Gene annotation: `ref/mm10.gtf`
   - ENCODE blacklist: `ref/mm10-blacklist.v2.bed`

Recommended input naming pattern:
- `SAMPLE_1.fastq.gz`
- `SAMPLE_2.fastq.gz`

This pattern is required by `fromFilePairs` in `main.nf`.

## Workflow overview

1. `FASTQC`: run FastQC for each paired-end sample.
2. `TRIM_GALORE` (optional): trim adapters and low-quality bases.
3. `ALIGN_BOWTIE2`: align to mm10 with `--very-sensitive`.
4. `FILTER_BAM`: MAPQ filter, remove `chrM`, remove blacklist regions, remove PCR duplicates.
5. `MACS2_PEAKS`: call open-chromatin peaks (`BAMPE`, `nomodel`, `keep-dup auto`, `extsize 147`).
6. `ATAC_QC_METRICS`: compute fragment distribution, FRiP, Tn5 shifted insertions, TSS profile/heatmap.
7. `MULTIQC`: aggregate FastQC/alignment/Picard/MACS2/trimming logs.

## Quick start

```bash
nextflow run main.nf -profile conda \
  --reads "data/raw_fastq/*_{1,2}.fastq.gz" \
  --bowtie2_index ref/bowtie2_index/mm10 \
  --blacklist ref/mm10-blacklist.v2.bed \
  --gtf ref/mm10.gtf \
  --macs2_gsize mm
```

Run with Singularity profile if needed:

```bash
nextflow run main.nf -profile singularity \
  --reads "data/raw_fastq/*_{1,2}.fastq.gz" \
  --bowtie2_index ref/bowtie2_index/mm10 \
  --blacklist ref/mm10-blacklist.v2.bed \
  --gtf ref/mm10.gtf
```

Resume an interrupted run:

```bash
nextflow run main.nf -profile conda -resume
```

## Main parameters

| Parameter | Default | Description |
|---|---|---|
| `--reads` | `data/raw_fastq/*_{1,2}.fastq.gz` | Input paired FASTQ glob |
| `--outdir` | `results` | Output root directory |
| `--bowtie2_index` | `ref/bowtie2_index/mm10` | Bowtie2 index prefix |
| `--gtf` | `ref/mm10.gtf` | Gene annotation for TSS-based QC |
| `--blacklist` | `ref/mm10-blacklist.v2.bed` | ENCODE blacklist regions |
| `--do_trim` | `true` | Enable/disable Trim Galore |
| `--min_mapq` | `30` | MAPQ threshold in BAM filtering |
| `--mito_name` | `chrM` | Mitochondrial chromosome name |
| `--macs2_gsize` | `mm` | MACS2 genome size |
| `--macs2_pval` | `0.01` | MACS2 p-value cutoff |
| `--macs2_extsize` | `147` | MACS2 extension size |
| `--tss_window` | `2000` | TSS flanking window size |
| `--nfr_max` | `100` | NFR upper fragment length cutoff |
| `--nbr_min` / `--nbr_max` | `180` / `247` | NBR fragment length range |

## Implemented ATAC-seq logic

- Raw read QC: FastQC + MultiQC
- Optional adapter trimming: Trim Galore
- Alignment: Bowtie2 with `--very-sensitive`
- BAM processing: sort, MAPQ filtering, mtDNA removal, blacklist removal, deduplication
- Peak calling: MACS2 (`--nomodel --keep-dup auto --extsize 147`)
- Advanced QC:
  - Tn5 shift correction (+4/-5)
  - Fragment size distribution
  - NFR/NBR split and TSS heatmap/profile
  - FRiP score

## Output structure

- `results/01_QC/FastQC`: per-sample FastQC reports
- `results/01_QC/MultiQC`: aggregated MultiQC report
- `results/01_QC/Trimmed`: trimming reports/logs
- `results/02_Alignment`: cleaned final BAM/BAI and Picard metrics
- `results/03_Peaks`: `narrowPeak`, `summits`, and MACS2 report tables
- `results/04_ATAC_QC`: FRiP, fragment size table, Tn5-shifted insertions, TSS plots

## Downstream analysis (R Markdown)

Use `notebooks/ATACseq_Analysis.Rmd` for:
- QC interpretation
- differential accessibility analysis (DiffBind/EdgeR)
- peak annotation (ChIPseeker)
- motif enrichment (motifmatchr + JASPAR2020)
- reproduction of Figure 6a–6f

## Reproducibility tips

- Keep `environment.yml` version-pinned.
- Use `-resume` to avoid recomputing completed tasks.
- Archive Nextflow run metadata (`timeline`, `report`, `trace`, `dag`) if enabled in your local profile.

## Notes

- Large data and intermediate files are excluded by `.gitignore`.
- Output folders are created under `results/01_QC`, `results/02_Alignment`, `results/03_Peaks`, and `results/04_ATAC_QC`.
- Notebook deliverables are in `notebooks/ATACseq_Analysis.Rmd`.
