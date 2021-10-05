#! /usr/bin/env python
# -*- coding utf-8 -*-

################################## INDEX REFS ##################################

rule star_index_gdc38_gencode38:
    conda: "envs/star.yaml"
    output:
        directory(config['indexes']['star'])
    input:
        genome = config['sequences']['genome'],
        annotation_gtf_gz = 'refs/downloads/gencode.v38.annotation.gtf.gz'
    params:
        sjdbOverhang = 74
    threads: workflow.cores
    shell:
        """
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)

        pigz -dc {input.gtf_gz} > $tdir/gencode.v38.annotation.gtf

        STAR\
            --runThreadN {threads}\
            --runMode genomeGenerate\
            --genomeDir {output}\
            --outFileNamePrefix {output}\
            --genomeFastaFiles {input.genome}\
            --sjdbGTFfile $tdir/gencode.v38.annotation.gtf\
            --sjdbOverhang {params.sjdbOverhang}
        """
