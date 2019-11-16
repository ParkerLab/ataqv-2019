#!/usr/bin/env python
# coding: utf-8

# In[55]:

import seaborn as sns
import os
import glob
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import re
import Ataqv
import itertools

# find the Su et al samples
SU_ET_AL_PROJECT = 'PRJNA323617.PRJNA341508'
sample_info = pd.read_csv(sys.argv[1], delimiter='\t')
sample_info = sample_info[sample_info.project==SU_ET_AL_PROJECT]
# some libraries have multiple runs...

TITLE_RE = 'GSM\d+: (.*)_rep(\d+) .*'
sample_info['condition'] = sample_info.title.map(lambda x: re.match(TITLE_RE, x).group(1))
sample_info['replicate'] = sample_info.title.map(lambda x: re.match(TITLE_RE, x).group(2))
# sample_info

flowcells = glob.glob('*.flowcells.txt')
def read_flowcell_file(f):
    tmp = None
    with open(f, 'r') as x:
        tmp = x.read().rstrip()
    return tmp


flowcells = {os.path.basename(f).replace('.flowcells.txt', ''): read_flowcell_file(f) for f in flowcells}
sample_info['flowcell'] = sample_info.run.map(lambda x: flowcells.get(x))
#sample_info


ataqv = list(itertools.chain(*[Ataqv.load_ataqv_json(f) for f in glob.glob('*.json.gz')]))

for a in ataqv:
    a.name = re.match('.*(SRR\d+)_.*', a.name).group(1)

metrics = [a for a in ataqv if a.name in list(sample_info.run)]
metrics = pd.DataFrame({metric: [getattr(a, metric) for a in ataqv] for metric in ['median_fragment_length', 'name']}).set_index('name')
#metrics.head()


sample_info['median_fragment_length'] = sample_info.run.map(lambda x: metrics.loc[x, 'median_fragment_length'])
p = sns.relplot(x='experiment', y='median_fragment_length', hue='condition', data=sample_info)
p = sns.relplot(x='flowcell', y='median_fragment_length', hue='condition', data=sample_info)
p = sns.relplot(x='experiment', y='median_fragment_length', hue='flowcell', data=sample_info)
p.fig.savefig('su_et_al_median_fragment_length_vs_batch_all_libraries.pdf')

keep_flowcells = sample_info[sample_info.condition.isin(['E0', 'E1'])].flowcell.unique()
p = sns.relplot(x='condition', y='median_fragment_length', hue='sequencing run/flowcell', data=sample_info[sample_info.flowcell.isin(keep_flowcells)].rename(columns={'flowcell': 'sequencing run/flowcell'}), alpha=0.3)
p.set_xlabels('Experimental condition')
p.set_ylabels('Median fragment length (bp)')
p.fig.savefig('su_et_al_median_fragment_length_vs_batch.pdf')
