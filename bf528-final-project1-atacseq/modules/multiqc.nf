process MULTIQC {
    tag "MULTIQC"
    publishDir "${params.outdir}/01_QC/MultiQC", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    mkdir -p multiqc_input
    for f in ${qc_files}; do
      cp -r \"$f\" multiqc_input/
    done

    multiqc multiqc_input -o .
    """
}
