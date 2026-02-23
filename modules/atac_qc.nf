process ATAC_QC_METRICS {
		tag "ATAC_QC on ${sample_id}"
		publishDir "${params.outdir}/04_ATAC_QC", mode: 'copy'

		input:
		tuple val(sample_id), path(bam), path(bai)
		tuple val(peak_sample_id), path(peaks)
		path gtf

		output:
		path "${sample_id}.fragment_size.tsv", emit: fragment_size
		path "${sample_id}.frip.txt", emit: frip
		path "${sample_id}.tn5_shifted.bed.gz", emit: tn5_shifted_bed
		path "${sample_id}.tss_enrichment.png", emit: tss_enrichment_plot
		path "${sample_id}.nfr_nbr_heatmap.png", emit: nfr_nbr_heatmap
		path "${sample_id}.tss_qc.log", emit: tss_log

		when:
		sample_id == peak_sample_id

		script:
		"""
		samtools stats ${bam} | awk 'BEGIN{OFS="\t"} /^IS/ {print \$2,\$3}' > ${sample_id}.fragment_size.tsv

		total_reads=\$(samtools view -@ ${task.cpus} -c -F 260 ${bam})
		reads_in_peaks=\$(bedtools intersect -u -abam ${bam} -b ${peaks} | samtools view -@ ${task.cpus} -c)
		awk -v t=\${total_reads} -v p=\${reads_in_peaks} 'BEGIN{frip=(t>0)?p/t:0; printf("sample_id\ttotal_reads\treads_in_peaks\tFRiP\\n%s\\t%d\\t%d\\t%.6f\\n", "${sample_id}", t, p, frip)}' > ${sample_id}.frip.txt

		bedtools bamtobed -i ${bam} \
			| awk 'BEGIN{OFS="\t"} {if(\$6=="+"){s=\$2+4; e=s+1}else{s=\$3-5; e=s+1} if(s<0){s=0; e=1} print \$1,s,e,\$4,\$5,\$6}' \
			| sort -k1,1 -k2,2n \
			| bgzip > ${sample_id}.tn5_shifted.bed.gz

		awk 'BEGIN{OFS="\t"} \$3=="gene" {if(\$7=="+"){tss=\$4-1}else{tss=\$5-1}; if(tss<0)tss=0; print \$1,tss,tss+1,"gene","0",\$7}' ${gtf} > ${sample_id}.tss.bed

		samtools view -h ${bam} \
			| awk 'BEGIN{OFS="\t"} /^@/ {print; next} {len=(\$9<0)?-\$9:\$9; if(len>0 && len<${params.nfr_max}) print}' \
			| samtools view -@ ${task.cpus} -bS -o ${sample_id}.nfr.bam -

		samtools view -h ${bam} \
			| awk 'BEGIN{OFS="\t"} /^@/ {print; next} {len=(\$9<0)?-\$9:\$9; if(len>=${params.nbr_min} && len<=${params.nbr_max}) print}' \
			| samtools view -@ ${task.cpus} -bS -o ${sample_id}.nbr.bam -

		samtools sort -@ ${task.cpus} -o ${sample_id}.nfr.sorted.bam ${sample_id}.nfr.bam
		samtools sort -@ ${task.cpus} -o ${sample_id}.nbr.sorted.bam ${sample_id}.nbr.bam
		samtools index -@ ${task.cpus} ${sample_id}.nfr.sorted.bam
		samtools index -@ ${task.cpus} ${sample_id}.nbr.sorted.bam

		bamCoverage -b ${sample_id}.nfr.sorted.bam -o ${sample_id}.nfr.bw --binSize ${params.bw_bin_size} --normalizeUsing CPM -p ${task.cpus}
		bamCoverage -b ${sample_id}.nbr.sorted.bam -o ${sample_id}.nbr.bw --binSize ${params.bw_bin_size} --normalizeUsing CPM -p ${task.cpus}

		computeMatrix reference-point \
			--referencePoint TSS \
			-b ${params.tss_window} -a ${params.tss_window} \
			-R ${sample_id}.tss.bed \
			-S ${sample_id}.nfr.bw ${sample_id}.nbr.bw \
			--missingDataAsZero \
			-o ${sample_id}.tss.matrix.gz \
			2> ${sample_id}.tss_qc.log

		plotProfile -m ${sample_id}.tss.matrix.gz \
			--samplesLabel NFR NBR \
			-out ${sample_id}.tss_enrichment.png

		plotHeatmap -m ${sample_id}.tss.matrix.gz \
			--samplesLabel NFR NBR \
			-out ${sample_id}.nfr_nbr_heatmap.png
		"""
}
