#! /usr/bin/env python
# -*- coding utf-8 -*-

################################### TELESCOPE ###################################

rule telescope:
    conda: "../envs/telescope.yaml"
    output:
        "results/telescope/{s}/{s}_telescope.report.tsv",
        "results/telescope/{s}/{s}_telescope.updated.bam",
        "results/telescope/{s}/{s}_telescope.other.bam"
    input:
        bam = "results/star_alignment/{s}/{s}_GDC38.Aligned.out.bam",
        annotation = rules.telescope_annotation.output
    benchmark: "benchmarks/telescope/{s}_telescope.tsv"
    log:
        "results/telescope/{s}/telescope.log"
    threads: config['telescope_threads']
    params:
        tmpdir = config['local_tmp']
    shell:
        """
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.s}.XXXXXX)
        telescope bulk assign\
         --exp_tag inform\
         --theta_prior 200000\
         --max_iter 200\
         --updated_sam\
         --outdir $tdir\
         {input[0]}\
         {input[1]}\
         2>&1 | tee {log[0]}
        mv $tdir/inform-TE_counts.tsv {output[0]}
        mv $tdir/inform-updated.bam {output[1]}
        mv $tdir/inform-other.bam {output[2]}
        chmod 660 {output[1]}
        rm -rf $tdir
        """

rule sample_complete:
    input:
        rules.telescope.output
    output:
        touch("results/completed/{s}_completed.txt")
