#! /usr/bin/env python
# -*- coding utf-8 -*-

# Add quant mode for gene TE_counts
# Make sure that star runs with stringtie
################################ STAR ALIGNMENT ################################

rule star_alignment:
    conda:
        "../envs/star.yaml"
    input:
        R1 = "samples/{s}/original_R1.fastq",
        R2 = "samples/{s}/original_R2.fastq",
        genome = config['indexes']['star']
    output:
        aligned_bam = "results/{s}/{s}_GDC38.Aligned.out.bam",
        sorted_bam = "results/{s}/{s}_GDC38.Aligned.sortedByCoord.out.bam"
    params:
        out_prefix="results/{s}/{s}_GDC38."
    threads: 8
    resources:
        mem_mb=70000, disk_mb=20000
    benchmark: "benchmarks/star_alignment/{s}_star_alignment.tsv"
    shell:
        '''
        STAR\
            --runThreadN {threads}\
            --genomeDir {input.genome}\
            --readFilesIn {input.R1} {input.R2}\
            --outSAMattributes NH HI NM MD AS XS\
            --outSAMtype BAM Unsorted SortedByCoordinate\
            --outFileNamePrefix {params.out_prefix}\
            --quantMode TranscriptomeSAM GeneCounts\
            --outSAMstrandField intronMotif\
            --limitBAMsortRAM 53679965568
        '''
