process FASTQC {
    tag "FASTQC on ${sample_id}"
    publishDir "${params.outdir}/01_QC/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
    path "*_fastqc.*", emit: qc_results

    script:
    """
    fastqc -t ${task.cpus} -q ${reads[0]} ${reads[1]}
    """
}
