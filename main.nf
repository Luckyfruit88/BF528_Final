nextflow.enable.dsl=2

params.accession_list = params.accession_list ?: "SRR_Acc_List.txt"
params.bowtie2_index  = params.bowtie2_index  ?: "path/to/mm10/index"
params.blacklist      = params.blacklist      ?: "path/to/blacklist.bed"
params.annotation_gtf = params.annotation_gtf ?: "path/to/annotation.gtf"
params.rna_counts     = params.rna_counts     ?: []
params.outdir         = params.outdir         ?: "results"

process DOWNLOAD_AND_DUMP {
    tag { accession }
    publishDir "${params.outdir}/raw", mode: 'copy'

    input:
    val accession

    output:
    path "${accession}.fastq.gz", emit: reads

    script:
    """
    prefetch ${accession}
    fastq-dump --split-files --gzip ${accession}/${accession}.sra
    mv ${accession}_1.fastq.gz ${accession}.fastq.gz
    """
}

process QC_AND_TRIM {
    tag { sample_id }
    publishDir "${params.outdir}/qc", mode: 'copy'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "${sample_id}_trimmed.fq.gz", emit: trimmed
    path "${sample_id}_fastqc.zip", emit: fastqc_zip
    path "${sample_id}_fastqc.html", emit: fastqc_html

    script:
    def sample_id = reads.baseName.replace('.fastq','').replace('.fq','')
    """
    fastqc ${reads}
    trim_galore --nextera --gzip ${reads}
    mv ${sample_id}_trimmed.fq.gz ${sample_id}_trimmed.fq.gz
    mv ${sample_id}_fastqc.zip ${sample_id}_fastqc.zip
    mv ${sample_id}_fastqc.html ${sample_id}_fastqc.html
    """
}

process ALIGN_BOWTIE2 {
    tag { sample_id }
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    path trimmed

    output:
    path "${sample_id}.sorted.bam", emit: bam
    path "${sample_id}.sorted.bam.bai", emit: bai

    script:
    def sample_id = trimmed.baseName.replace('_trimmed','')
    """
    bowtie2 --very-sensitive -x ${params.bowtie2_index} -U ${trimmed} | \
        samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process FILTER_ALIGNMENTS {
    tag { sample_id }
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    path bam

    output:
    path "${sample_id}_final.bam", emit: filtered_bam
    path "${sample_id}_final.bam.bai", emit: filtered_bai

    script:
    def sample_id = bam.baseName.replace('.sorted','')
    """
    samtools view -h -q 30 -F 4 ${bam} | \
        grep -v "chrM" | \
        bedtools intersect -v -a stdin -b ${params.blacklist} | \
        samtools sort -o ${sample_id}_final.bam
    samtools index ${sample_id}_final.bam
    """
}

process CALL_PEAKS {
    tag { sample_id }
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    path filtered_bam

    output:
    path "${sample_id}_peaks.narrowPeak", emit: peaks

    script:
    def sample_id = filtered_bam.baseName.replace('_final','')
    """
    macs2 callpeak -t ${filtered_bam} -f BAM -g mm --nomodel --extsize 147 --keep-dup auto -q 0.01 -n ${sample_id}_peaks
    """
}

process DIFFERENTIAL_ACCESSIBILITY {
    tag "diffbind"
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path peaks
    path final_bams

    output:
    tuple path('DARs_cDC1.csv'), path('DARs_cDC2.csv'), emit: dars_pair

    script:
    def peak_env = peaks.join(';')
    def bam_env  = final_bams.join(';')
    """
    export PEAK_FILES="${peak_env}"
    export BAM_FILES="${bam_env}"

    Rscript - <<'RSCRIPT'
    peaks <- strsplit(Sys.getenv('PEAK_FILES'), ';')[[1]]
    bams  <- strsplit(Sys.getenv('BAM_FILES'), ';')[[1]]

    sample_sheet <- data.frame(
      SampleID = sub('\\\.bam$', '', basename(bams)),
      Condition = ifelse(grepl('KO', bams, ignore.case = TRUE), 'KO', 'WT'),
      Replicate = seq_along(bams),
      bamReads = bams,
      Peaks = peaks,
      stringsAsFactors = FALSE
    )
    write.csv(sample_sheet, 'diffbind_samples.csv', row.names = FALSE)

    suppressPackageStartupMessages({
      library(DiffBind)
      library(edgeR)
    })

    db_obj <- dba(sampleSheet = 'diffbind_samples.csv')
    db_obj <- dba.count(db_obj, bUseSummarizeOverlaps = TRUE)
    db_obj <- dba.normalize(db_obj, method = DBA_EDGER, normalize = DBA_NORM_TMM)
    db_obj <- dba.contrast(db_obj, categories = DBA_CONDITION)
    db_obj <- dba.analyze(db_obj, method = DBA_EDGER)

    res <- dba.report(db_obj)
    res_df <- as.data.frame(res)

    cdc1 <- res_df[grepl('cDC1', res_df$Conc, ignore.case = TRUE) | grepl('cDC1', res_df$Gene.Name, ignore.case = TRUE), ]
    cdc2 <- res_df[grepl('cDC2', res_df$Conc, ignore.case = TRUE) | grepl('cDC2', res_df$Gene.Name, ignore.case = TRUE), ]

    write.csv(cdc1, 'DARs_cDC1.csv', row.names = FALSE)
    write.csv(cdc2, 'DARs_cDC2.csv', row.names = FALSE)
    RSCRIPT
    """
}

process INTEGRATE_MULTIOMICS {
    tag "integration"
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    tuple path(dars_cdc1), path(dars_cdc2)
    val rna_counts

    output:
    path 'integration_outputs', emit: integration_dir

    script:
    def rna_args = (rna_counts instanceof List) ? rna_counts.join(' ') : rna_counts
    """
    Rscript - <<'RSCRIPT' ${rna_args}
    args <- commandArgs(trailingOnly = TRUE)
    dars <- c('${dars_cdc1}', '${dars_cdc2}')
    rna_files <- args

    suppressPackageStartupMessages({
      library(DESeq2)
      library(ChIPseeker)
      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      library(org.Mm.eg.db)
      library(clusterProfiler)
      library(motifmatchr)
      library(JASPAR2020)
      library(TFBSTools)
      library(ggplot2)
      library(data.table)
    })

    dir.create('integration_outputs', showWarnings = FALSE)

    dars_list <- lapply(dars, data.table::fread)
    names(dars_list) <- c('cDC1', 'cDC2')

    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peak_gr <- lapply(dars_list, function(df) GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE))
    annotated <- lapply(peak_gr, function(gr) annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = 'org.Mm.eg.db'))

    count_list <- lapply(rna_files, function(f) data.table::fread(f))
    names(count_list) <- paste0('rna', seq_along(count_list))

    deseq_results <- lapply(count_list, function(dt){
      gene_ids <- dt[[1]]
      counts <- as.matrix(dt[,-1])
      condition <- factor(c(rep('WT', ncol(counts)/2), rep('KO', ncol(counts)/2)))
      coldata <- data.frame(condition = condition)
      dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
      dds <- DESeq(dds)
      res <- results(dds)
      res$gene <- gene_ids
      as.data.frame(res)
    })

    for (celltype in names(annotated)) {
      anno_df <- as.data.frame(annotated[[celltype]])
      integration_tables <- lapply(deseq_results, function(res) merge(anno_df, res, by.x = 'geneId', by.y = 'gene', all = FALSE))
      combined <- integration_tables[[1]]
      if (length(integration_tables) > 1) {
        for (i in 2:length(integration_tables)) {
          combined <- merge(combined, integration_tables[[i]], by = intersect(names(combined), names(integration_tables[[i]])), all = FALSE)
        }
      }
      write.csv(combined, file.path('integration_outputs', paste0('Integrated_', celltype, '.csv')), row.names = FALSE)

      p <- ggplot(combined, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
        geom_point(alpha = 0.5) +
        labs(title = paste0('Accessibility vs Expression: ', celltype), x = 'Accessibility log2FC', y = 'Expression log2FC') +
        theme_minimal()
      ggsave(file.path('integration_outputs', paste0('Scatter_', celltype, '.pdf')), p, width = 6, height = 5)

      go_res <- enrichGO(gene = unique(combined$geneId), OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = 'BP')
      pdf(file.path('integration_outputs', paste0('GO_', celltype, '.pdf')))
      print(dotplot(go_res, showCategory = 20))
      dev.off()

      motifs <- getMatrixSet(x = JASPAR2020, opts = list(species = 10090, all_versions = FALSE))
      motif_hits <- matchMotifs(motifs, peak_gr[[celltype]])
      motif_counts <- colSums(SummarizedExperiment::assay(motif_hits))
      motif_df <- data.frame(motif = names(motif_counts), count = motif_counts)
      motif_df <- motif_df[order(-motif_df$count), ][1:min(20, nrow(motif_df)), ]
      pdf(file.path('integration_outputs', paste0('Motifs_', celltype, '.pdf')))
      barplot(motif_df$count, names.arg = motif_df$motif, las = 2, main = paste0('Top motifs: ', celltype))
      dev.off()
    }
    RSCRIPT
    """
}

process TRACKS_AND_QC {
    tag { sample_id }
    publishDir "${params.outdir}/tracks", mode: 'copy'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path filtered_bam

    output:
    path "${sample_id}.bigWig", emit: bigwig
    path "${sample_id}_TSS_heatmap.pdf", emit: heatmap

    script:
    def sample_id = filtered_bam.baseName.replace('_final','')
    """
    bamCoverage -b ${filtered_bam} -o ${sample_id}.bigWig --normalizeUsing CPM --binSize 10 --smoothLength 30 --ignoreForNormalization chrX chrY chrM
    computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -S ${sample_id}.bigWig -R ${params.annotation_gtf} -o ${sample_id}_tss_matrix.gz
    plotHeatmap -m ${sample_id}_tss_matrix.gz -out ${sample_id}_TSS_heatmap.pdf
    """
}

workflow {
    Channel.fromPath(params.accession_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .set { accessions }

    downloaded = DOWNLOAD_AND_DUMP(accessions)
    qc_trimmed = QC_AND_TRIM(downloaded)
    aligned = ALIGN_BOWTIE2(qc_trimmed)
    filtered = FILTER_ALIGNMENTS(aligned)

    def filt_split = filtered.out.filtered_bam.into { bam_for_peaks; filt_for_diff; filt_for_tracks }
    peaks_called = CALL_PEAKS(bam_for_peaks)
    def peaks_split = peaks_called.out.peaks.into { peaks_for_diff; peaks_for_tracks }

    def diff_peaks = peaks_for_diff.collect()
    def diff_bams  = filt_for_diff.collect()
    differential_results = DIFFERENTIAL_ACCESSIBILITY(diff_peaks, diff_bams)

    def rna_counts_ch = Channel.value(params.rna_counts)
    INTEGRATE_MULTIOMICS(differential_results.out.dars_pair, rna_counts_ch)

    TRACKS_AND_QC(filt_for_tracks)
}
