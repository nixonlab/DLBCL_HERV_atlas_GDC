#! /usr/bin/env python
# -*- coding: utf-8 -*-

################################# DOWNLOAD BAMS #################################


rule gdc_download:
    """ Download BAM from Genomic Data Commons (GDC)
    """
    input:
        token = config['gdc_token_file']
    output:
        temp("samples/{s}/original.bam")
    params:
        uuid = lambda wc: PILOT[wc.s]['file_uuid'],
        md5sum = lambda wc: PILOT[wc.s]['md5sum']
    shell:
        '''
mkdir -p $(dirname {output[0]})
curl\
 -H "X-Auth-Token: $(<{input[0]})"\
 https://api.gdc.cancer.gov/data/{params.uuid}\
 > {output[0]}
echo {params.md5sum} {output[0]} | md5sum -c -
chmod 600 {output[0]}
        '''

rule download_complete:
    input:
        expand("samples/{s}/original.bam", s = PILOT.keys())
