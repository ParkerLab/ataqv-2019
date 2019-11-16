#!/usr/bin/env python
# coding: utf-8

# In[8]:

import os
import sys
import argparse
import pybedtools as bt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
import re
import logging
from adjustText import adjust_text
import gzip
import csv
import itertools
import scipy.stats as ss
from scipy.stats.stats import pearsonr

plt.ioff()

parser = argparse.ArgumentParser(description='Calculate TSS enrichment from a bam file.', add_help = True)
parser.add_argument('--sample-info', dest='sample_info', type = str,  help = 'BAM file, sorted by read name.')
parser.add_argument('--ataqv-metrics', dest = 'ataqv_metrics', help = 'Ataqv metrics')
parser.add_argument('--inverse-rank-normalization', dest = 'irn', action='store_true', help = 'Perform inverse rank normalization on the metrics.')

args = parser.parse_args()


def inverse_rank_norm(values):
    """Perform rank-based inverse normal transform to a numeric Pandas series"""
    #values = pd.Series([5, 7, 2, 1, 1])
    quantiles = (values.rank()-0.5)/(len(values))
    return ss.norm.ppf(quantiles)


# In[10]:

sample_info = pd.read_csv(args.sample_info, delimiter='\t')
sample_info = sample_info.loc[~sample_info.is_single_cell,['experiment', 'organism', 'project']].drop_duplicates().set_index('experiment')


# In[34]:

ataqv = pd.read_csv(args.ataqv_metrics, delimiter='\t', header=None)
ataqv.columns = ['experiment', 'metric', 'value']
ataqv.value = ataqv.value.map(lambda x: np.nan if x == 'None' else float(x))


# In[35]:

metrics = ataqv.pivot(columns='metric', values='value', index='experiment')


# In[38]:

#metrics.head()


# In[37]:

# drop cols with all 0s
metrics = metrics.loc[:,metrics.apply(lambda x: any(x>0), axis='index').values]


# In[39]:

metrics.fillna(0, inplace=True)
metrics.replace(-np.inf, 0, inplace=True)
#metrics.head()


# In[15]:

# inverse rank norm
if args.irn:
    for i in list(metrics.columns):
        df[[metrics]] = inverse_rank_norm(df[[metrics]])


# In[40]:

correlation_significance = pd.DataFrame(columns=metrics.columns, index=metrics.columns)
for i in list(metrics.columns):
    for j in list(metrics.columns):
        correlation_significance.loc[i,j] = pearsonr(metrics[i].values, metrics[j].values)[1]
#correlation_significance = correlation_significance.applymap(lambda x: '*' if x < 0.01 else '')
#correlation_significance.head()
correlation_significance.to_csv('all-correlation-p-values.txt', sep='\t')


# In[43]:

#plt.rcParams['figure.figsize'] = (8, 7)
#p = sns.clustermap(metrics.corr(), cmap='bwr', vmax=1, vmin=-1, cbar_kws={'label': 'Pearson corr.'})


# In[73]:

# heatmap versions of pairplots
plt.rcParams['figure.figsize'] = (8, 7)
p = sns.clustermap(metrics.corr(), cmap='bwr', vmax=1, vmin=-1, cbar_kws={'label': 'Pearson corr.'})
p.savefig(f'all-corr-heatmap.pdf')
p.fig.clear()

