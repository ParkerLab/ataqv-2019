
# coding: utf-8

# In[24]:

import os
import sys
import argparse
import pybedtools as bt
import numpy as np
import pandas as pd
import re
import logging
import gzip
import csv
from Ataqv import Ataqv, load_ataqv_json
import itertools
import glob

parser = argparse.ArgumentParser(description='Create covariate file for DESeq script (from sample info file and ataqv metric files).', add_help = True)
parser.add_argument('sample_info', type = str,  help = 'Sample info file.')
parser.add_argument('--metric-files', dest='metric_files', nargs='+', type = str, help = 'Ataqv json files.')
parser.add_argument('--include', dest='include', nargs='+', type = str, help = 'List of variable to include (e.g., library, tss_enrichment, replicate, ...). May be an ataqv metric or a column from the sample info file.')
parser.add_argument('--libraries', dest='libraries', nargs='+', type = str, help = 'List of libraries to keep.')

args = parser.parse_args()


sample_info = pd.read_csv(args.sample_info, delimiter='\t').set_index('library')
sample_info.index = sample_info.index.map(str)


# In[26]:

ataqv = list(itertools.chain(*[load_ataqv_json(f) for f in args.metric_files]))

ataqv_metrics = [i for i in args.include if i not in sample_info.columns] + ['library']
ataqv_metrics = list(set(ataqv_metrics))
ataqv_metrics
# verify that the metrics are present...
for i in ataqv_metrics:
    assert(hasattr(ataqv[0], i))

ataqv = pd.DataFrame({metric: [getattr(a, metric) for a in ataqv] for metric in ataqv_metrics}).set_index('library')


# In[37]:

df = sample_info.join(ataqv).reset_index()
df = df.loc[df.library.isin(args.libraries),args.include].drop_duplicates()
df.to_csv(sys.stdout, sep='\t', index = False, header = True)


# In[ ]:



