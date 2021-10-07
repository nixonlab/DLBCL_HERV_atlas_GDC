#! /usr/bin/env python
# -*- coding utf-8 -*-

################################ STAR ALIGNMENT ################################

rule star_alignment:
	conda: "envs/star.yaml"
	input:
		R1 = "samples/{sampid}/original_R1.fastq",
		R2 = "samples/{sampid}/original_R2.fastq",
	output:
		aligned_bam = "results/{sampid}_GDC38.Aligned.out.bam",
        sorted_bam = "results/{sampid}_GDC38.Aligned.sortedByCoord.out.bam"
	params:
		genome = config['indexes']['star'],
		out_prefix = "results/{sampid}_GDC38."
	threads: 18
	shell:
		'''
STAR --runThreadN {threads}\
--genomeDir {params.genome}\
--readFilesIn {input.R1},{input.R2}
--outSAMattributes All\
--outSAMtype BAM Unsorted SortedByCoordinate\
--outFileNamePrefix {params.out_prefix}
		'''
