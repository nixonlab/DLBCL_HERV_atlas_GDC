#! /usr/bin/env python
# -*- coding utf-8 -*-

################################## INDEX REFS ##################################

rule star_index:
    conda: "envs/star.yaml"
    output:
        directory(config['indexes']['star'])
    input:
        genome = config['sequences']['genome'],
        annotation_gtf_gz = config['annotations']['gencode']
    params:
        sjdbOverhang = 74
    threads: 8
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
