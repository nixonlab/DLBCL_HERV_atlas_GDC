This folder is meant to contain all resources necessary for running the workflow, for example reference sequences or databases.
Wherever feasible, they can also be downloaded programmatically via rules defined in the pipeline.

### Sample metadata

The script `gdcquery.py` performs a query through the GDC API files endpoint and retrieves
all RNA-seq BAM files belonging to a specific project.


###### [NCICCR-DLBCL sample matrix](resources/metadata/NCICCR-DLBCL.tsv)

The NCI CCR DLBCL project includes three types of alignments: genomic, transcriptomic, and
chimeric. ("Chimeric" is produced by STAR and contains chimeric alignments only). Since we
are going to be extracting and remapping reads, the genomic alignments are the ones we 
want. Here is how we got the sample matrix:

```bash
workflow/scripts/gdcquery.py NCICCR-DLBCL > resources/metadata/tmp.txt
head -n1 resources/metadata/tmp.txt > resources/metadata/NCICCR-DLBCL.tsv
grep 'genomic' resources/metadata/tmp.txt >> resources/metadata/NCICCR-DLBCL.tsv
rm resources/metadata/tmp.txt 
```

##### [TCGA-DLBC sample matrix](resources/metadata/TCGA-DLBC.tsv)

```bash
workflow/scripts/gdcquery.py TCGA-DLBC > resources/metadata/TCGA-DLBC.tsv
```


##### [Table S9 DLBCL Characteristics](resources/metadata/tableS9.tsv)

DLBCLs have been previously characterized by Schmitz et al. Genetics and Pathogenesis of Diffuse Large B-Cell Lymphoma. 
The excel file was downloaded from [https://www.nejm.org/doi/suppl/10.1056/NEJMoa1801445/suppl_file/nejmoa1801445_appendix_2.xlsx](https://www.nejm.org/doi/suppl/10.1056/NEJMoa1801445/suppl_file/nejmoa1801445_appendix_2.xlsx)
Table S9 was copied from the excel file into a tsv.

##### 

The script `format_samples.py` renames columns in `TCGA-DLBC.tsv` and `NCICCR-DLBCL.tsv`
and appends gene expression subgroup and genetic subtype information to each sample.
Pilot sets are randomly selected from the NCICCR project with up to 10 samples from each 
Gene Expression X Genetic Subtype combination (see below), for a total of 98 samples.
The resulting sample sheet is saved as [`config/samples.tsv`](./config/samples.tsv).

```bash
python workflow/scripts/format_samples.py
```

The number of samples for each Gene Expression X Genetic Subtype combination. Up to 10 
samples randomly selected from each combination selected for pilot. 

```bash
Gene Expression X Genetic Subtype counts (NCICCR only)
('ABC', 'Other')	138
('GCB', 'Other')	64
('Unclass', 'Other')	59
('GCB', 'EZB')	57
('ABC', 'MCD')	54
('ABC', 'BN2')	36
('Unclass', 'BN2')	35
('GCB', 'BN2')	16
('ABC', 'N1')	14
('Unclass', 'EZB')	4
('Unclass', 'MCD')	1
('Unclass', 'N1')	1
('ABC', 'EZB')	1
('GCB', 'MCD')	1
```