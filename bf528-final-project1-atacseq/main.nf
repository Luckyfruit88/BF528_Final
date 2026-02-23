#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC }        from './modules/fastqc.nf'
include { TRIM_GALORE }   from './modules/trim_galore.nf'
include { ALIGN_BOWTIE2 } from './modules/bowtie2.nf'
include { FILTER_BAM }    from './modules/filter_bam.nf'
include { MACS2_PEAKS }   from './modules/macs2.nf'
include { ATAC_QC_METRICS } from './modules/atac_qc.nf'
include { MULTIQC }       from './modules/multiqc.nf'

workflow {
    // Build paired-end input channel
    raw_pairs = Channel
        .fromFilePairs(params.reads, size: 2, checkIfExists: true)
        .ifEmpty { error "No FASTQ pairs were found with pattern: ${params.reads}" }

    // Resolve static reference inputs as value channels
    idx_ch = Channel.value(file(params.bowtie2_index))
    bl_ch  = Channel.value(file(params.blacklist))
    gtf_ch = Channel.value(file(params.gtf))

    FASTQC(raw_pairs)

    // Optional trimming branch, always producing a log channel
    if (params.do_trim) {
        TRIM_GALORE(raw_pairs)
        align_ch = TRIM_GALORE.out.trimmed_reads
        trim_logs_ch = TRIM_GALORE.out.trim_logs
    } else {
        align_ch = raw_pairs.map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
        trim_logs_ch = Channel.empty()
    }

    ALIGN_BOWTIE2(align_ch, idx_ch)

    FILTER_BAM(ALIGN_BOWTIE2.out.sorted_bam, bl_ch)

    MACS2_PEAKS(FILTER_BAM.out.final_bam)

    ATAC_QC_METRICS(
        FILTER_BAM.out.final_bam,
        MACS2_PEAKS.out.peaks,
        gtf_ch
    )

    qc_files = FASTQC.out.qc_results
        .mix(trim_logs_ch)
        .mix(ALIGN_BOWTIE2.out.align_log)
        .mix(FILTER_BAM.out.picard_metrics)
        .mix(MACS2_PEAKS.out.macs2_log)

    MULTIQC(qc_files.collect())
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
