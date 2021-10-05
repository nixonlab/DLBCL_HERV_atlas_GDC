#! /usr/bin/env python
# -*- coding utf-8 -*-

################################### TELESCOPE ###################################

rule telescope:
    conda: "envs/telescope"
    output:
    input:
        bam = "results/{sample}_GDC38.Aligned.out.bam",
        gtf = config['annotations']['gencode']
    params:
        tmpdir = config['local_tmp']
    log: "results/{sample}_telescope.log"
    shell:
    
