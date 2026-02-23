process MACS2_PEAKS {
    tag "MACS2 on ${sample_id}"
    publishDir "${params.outdir}/03_Peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: peaks
    path "${sample_id}_peaks.narrowPeak", emit: narrow_peak
    path "${sample_id}_peaks.xls", emit: peak_xls
    path "${sample_id}_summits.bed", emit: summits
    path "${sample_id}.macs2.log", emit: macs2_log

    script:
    """
    macs2 callpeak \
        -t ${bam} \
        -f BAMPE \
        -n ${sample_id} \
        -g ${params.macs2_gsize} \
        -p ${params.macs2_pval} \
        --nomodel \
        --keep-dup auto \
        --extsize ${params.macs2_extsize} \
        2> ${sample_id}.macs2.log
    """
}
