######################################################
## Analyze droplet single-cell RNA sequencing data
## with STARsolo
######################################################
## to do:
## input as function for cDNA and barcodes
## a scheme to name output directories
## see if I can integrate PEP here

rule starsolo_alignment:
	"""
	Align sequencing reads from a 10x V3 single-cell RNA-seq experiment using STARsolo
	"""
	input:
		R1 = "samples/{sampid}/original_R1.fastq",
		R2 = "samples/{sampid}/original_R2.fastq",
		genome = directory(config['indexes']['star']),
		whitelist = # FILL THIS IN
	output:
		"results/{sample}_GDC38.Aligned.out.bam",
        "results/{sample}_GDC38.Aligned.sortedByCoord.out.bam"
	params:
		out_prefix="results/{sample}_GDC38.",
		cb_start = config['cellbarcode_start'], # FILL THIS IN
		cb_length = config['cellbarcode_length'], # FILL THIS IN
		umi_start = config['umi_start'], # FILL THIS IN
		umi_length = config['umi_length'], # FILL THIS IN
		max_multimap = config['max_multimap'] # FILL THIS IN
	conda:
		"../envs/star.yaml"
	threads: 18
	shell:
		'''
		#--- STARsolo (turned on by --soloType CB_UMI_Simple)
		STAR\
			--runThreadN {threads}\
			--genomeDir {input.genome}\
			--readFilesIn {input.cDNA_L1},{input.cDNA_L2} {input.barcode_L1},{input.barcode_L2}\
			--readFilesCommand gunzip -c\
			--soloType CB_UMI_Simple\
			--soloCBwhitelist {input.whitelist}\
			--soloCBstart {params.cb_start}\
			--soloCBlen {params.cb_length}\
			--soloUMIstart {params.umi_start}\
			--soloUMIlen {params.umi_length}\
			--outFilterMultimapNmax {params.max_multimap}\
			--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
			--outSAMtype BAM Unsorted SortedByCoordinate\
			--outFileNamePrefix {params.out_prefix}
		'''
