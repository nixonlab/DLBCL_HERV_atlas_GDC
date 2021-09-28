#! /usr/bin/env python

from collections import Counter
from collections import OrderedDict
import random

# Reformat TCGA-DLBC
tmp1 = []
tsv = (l.strip('\n').split('\t') for l in open('resources/metadata/TCGA-DLBC.tsv', 'r'))
qhead = next(tsv)
for row in tsv:
    tmp1.append(dict(zip(qhead, row)))

remap = OrderedDict({
    'sample': 'cases.0.samples.0.submitter_id',
    'project': 'cases.0.project.project_id',
    'case': 'cases.0.submitter_id',
    'file': 'file_name',
    'md5sum': 'md5sum',
    'case_uuid': 'cases.0.case_id',
    'sample_uuid': 'cases.0.samples.0.sample_id',
    'file_uuid': 'file_id',
})

tmp1a = [{k:row[v] for k,v in remap.items()} for row in tmp1]

# Reformat NCICCR-DLBCL
tmp2 = []
tsv = (l.strip('\n').split('\t') for l in open('resources/metadata/NCICCR-DLBCL.tsv', 'r'))
qhead = next(tsv)
for row in tsv:
    tmp2.append(dict(zip(qhead, row)))

tmp2a = [{k:row[v] for k,v in remap.items()} for row in tmp2]

# take off the -sample from sample name
for row in tmp2a:
    row['sample'] = row['sample'].replace('-sample', '')


# Create combined samples
samples = tmp1a + tmp2a

# Load supplemental table S9 from Schmitz et al. (doi:10.1056/NEJMoa1801445)
tableS9 = {}
tsv = (l.strip('\n').split('\t') for l in open('resources/metadata/tableS9.tsv', 'r'))
s9header = next(tsv)
for row in tsv:
    tableS9[row[0]] = dict(zip(s9header, row))

tableS5 = {}
tsv = (l.strip('\n').split('\t') for l in open('resources/metadata/tableS5.tsv', 'r'))
s5header = next(tsv)
for row in tsv:
    tableS5[row[0]] = dict(zip(s5header, row))


# Add Gene expression group and Genetic group to 
for row in samples:
    if row['project'] == 'NCICCR-DLBCL':
        assert row['sample'] in tableS9, "row %s not in table" % row['sample']
        row['gex_group'] = tableS9[row['sample']]['Gene Expression Subgroup']
        row['gen_group'] = tableS9[row['sample']]['Genetic Subtype']
    elif row['project'] == 'TCGA-DLBC':
        if row['case'] in tableS9:
            row['gex_group'] = tableS9[row['case']]['Gene Expression Subgroup']
            row['gen_group'] = tableS9[row['case']]['Genetic Subtype']
        else:
            print('not in S9: %s' % row['case'])
            row['gex_group'] = ''
            row['gen_group'] = ''


# Add Gene expression group and Genetic group to 
for row in samples:
    if row['sample'] in tableS5:
        row['scCOO_group'] = tableS5[row['sample']]['sc-COO Group']
    else:
        row['scCOO_group'] = ''
    print('%s\t%s' % (row['sample'], row['scCOO_group']))
    
Counter(row)


# select pilot (only from NCICCR)
print('Gene Expression Subgroup counts (NCICCR only)')
print('\n'.join('%s\t%d' % t for t in Counter(r['gex_group'] for r in tmp2a).most_common()))
print('Genetic Subtype counts (NCICCR only)')
print('\n'.join('%s\t%d' % t for t in Counter(r['gen_group'] for r in tmp2a).most_common()))
print('Gene Expression X Genetic Subtype counts (NCICCR only)')
print('\n'.join('%s\t%d' % t for t in Counter((r['gex_group'],r['gen_group']) for r in tmp2a).most_common()))
print('scCOO Group')
print('\n'.join('%s\t%d' % t for t in Counter(r['scCOO_group'] for r in tmp2a).most_common()))

random.seed(12345)
pilot = set()
ssize = 10
tot = 0
# combs = [_[0] for _ in Counter((r['gex_group'],r['gen_group']) for r in tmp2a).most_common()]
groups = ['Group I', 'Group II', 'Group III', 'Group IV', 'Group V']
for g in groups:
    sel = [r for r in samples if r['scCOO_group'] == g]
    if len(sel) > ssize:
        selid = set(r['sample'] for r in random.sample(sel, 10))
    else:
        selid = set(r['sample'] for r in sel)
    print('%s\t%d\t%d' % (g, len(sel), len(selid)))
    tot += len(selid)
    pilot |= selid
    
    
# for t in combs:
#     sel = [r for r in tmp2a if r['gex_group'] == t[0] and r['gen_group'] == t[1]]
#     if len(sel) > ssize:
#         selid = set(r['sample'] for r in random.sample(sel, 10))
#     else:
#         selid = set(r['sample'] for r in sel)
#     print('%s\t%d\t%d' % (t, len(sel), len(selid)))
#     tot += len(selid)
#     pilot |= selid

assert len(pilot) == tot
print('%d samples in pilot' % len(pilot))
print(len(pilot))

for row in samples:
    if row['sample'] in pilot:
        row['pilot'] = 'True'
    else:
        row['pilot'] = 'False'

newhead = list(remap.keys()) + ['gex_group', 'gen_group', 'scCOO_group', 'pilot']

with open('config/samples.tsv', 'w') as outh:
    print('\t'.join(newhead), file=outh)
    for row in samples:
        print('\t'.join(row[f] for f in newhead), file=outh)
  

