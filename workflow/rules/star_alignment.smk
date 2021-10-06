#! /usr/bin/env python
# -*- coding utf-8 -*-

################################ STAR ALIGNMENT ################################

rule star_alignment:
	input:
		R1 = "samples/{sampid}/original_R1.fastq",
		R2 = "samples/{sampid}/original_R2.fastq",
		genome = directory(config['indexes']['star'])
	output:
		aligned_bam = "results/{sampid}_GDC38.Aligned.out.bam",
        sorted_bam = "results/{sampid}_GDC38.Aligned.sortedByCoord.out.bam"
	params:
		out_prefix = "results/{sampid}_GDC38."
	conda: "envs/star.yaml"
	threads: 18
	shell:
		"""
		STAR\
			--runThreadN {threads}\
			--genomeDir {input.genome}\
			--readFilesIn {input.R1},{input.R2}
			--outSAMattributes All\
			--outSAMtype BAM Unsorted SortedByCoordinate\
			--outFileNamePrefix {params.out_prefix}
		"""
