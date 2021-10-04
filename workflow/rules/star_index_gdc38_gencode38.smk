rule star_index_gdc38_gencode38:
    input:
        genome = config['sequences']['genome'],
        annotation_gtf_gz = 'refs/downloads/gencode.v38.annotation.gtf.gz'
    output:
        directory(config['indexes']['star'])
    params:
        sjdbOverhang = 74
    conda: "../envs/star.yaml"
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
