process FILTER_BAM {
    tag "FILTER on ${sample_id}"
    publishDir "${params.outdir}/02_Alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai)
    path blacklist_bed

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: final_bam
    path "${sample_id}.picard.metrics.txt", emit: picard_metrics

    script:
    """
    samtools view -@ ${task.cpus} -h -q ${params.min_mapq} ${sorted_bam} | \
      awk 'BEGIN{OFS="\t"} /^@/ {print; next} \$3 != "${params.mito_name}" {print}' | \
      samtools view -@ ${task.cpus} -bS -o ${sample_id}.nomt.bam -

    bedtools intersect -v -abam ${sample_id}.nomt.bam -b ${blacklist_bed} > ${sample_id}.noblacklist.bam
    samtools sort -@ ${task.cpus} -o ${sample_id}.noblacklist.sorted.bam ${sample_id}.noblacklist.bam

    picard MarkDuplicates \
        I=${sample_id}.noblacklist.sorted.bam \
        O=${sample_id}.final.bam \
        M=${sample_id}.picard.metrics.txt \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true

    samtools index -@ ${task.cpus} ${sample_id}.final.bam
    """
}
