process ALIGN_BOWTIE2 {
    tag "BOWTIE2 on ${sample_id}"
    publishDir "${params.outdir}/02_Alignment", mode: 'copy', pattern: "*.bowtie2.log"

    input:
    tuple val(sample_id), path(read1), path(read2)
    path bowtie2_index_prefix

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam
    path "${sample_id}.bowtie2.log", emit: align_log

    script:
    """
    bowtie2 -p ${task.cpus} \
        --very-sensitive \
        -X 2000 \
        --no-mixed \
        --no-discordant \
        -x ${bowtie2_index_prefix} \
        -1 ${read1} \
        -2 ${read2} \
        2> ${sample_id}.bowtie2.log \
        | samtools view -@ ${task.cpus} -bS - \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """
}
