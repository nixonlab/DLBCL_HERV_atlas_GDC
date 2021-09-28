#! /usr/bin/env python

from collections import Counter
from collections import OrderedDict
import pandas as pd

samples_tsv = 'config/samples.tsv'
read_samples_tsv = pd.read_csv(samples_tsv, sep='\t')
