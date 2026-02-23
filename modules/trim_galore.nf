process TRIM_GALORE {
    tag "TRIM on ${sample_id}"
    publishDir "${params.outdir}/01_QC/Trimmed", mode: 'copy',
        pattern: "*.txt"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: trim_logs

    script:
    """
    trim_galore --paired --cores ${task.cpus} --gzip --basename ${sample_id} ${reads[0]} ${reads[1]}
    """
}
