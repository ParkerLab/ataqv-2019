#!/usr/bin/env python
# coding: utf-8

import pybedtools as bt
import numpy as np
import pandas as pd
import sys
import os
from matplotlib import pyplot as plt
import seaborn as sns
import itertools
import re


SAMPLE_INFO = sys.argv[1]
MATRIX = sys.argv[2]
PREFIX = sys.argv[3]

sample_info = pd.read_csv(SAMPLE_INFO, delimiter='\t')[['tn5', 'replicate', 'library']]
sample_info.library = sample_info.library.map(str)
sample_info = sample_info.set_index('library')
sample_info['label'] = sample_info[['tn5', 'replicate']].apply(lambda x: '{} Tn5, {}'.format(*x), axis='columns')

def peak_file_to_library(f):
    return re.match('^(.*)_peaks.broadPeak.*', os.path.basename(f)).group(1)


df = pd.read_csv(MATRIX, delimiter='\t')
df['library_1'] = df.file_1.map(peak_file_to_library)
df['library_2'] = df.file_2.map(peak_file_to_library)
# df = df.drop(columns=['file_1', 'file_2'])
df = df.drop(['file_1', 'file_2'], axis=1)


LIBRARIES = df.library_1.append(df.library_2).unique()

def pair2jaccard(library_1, library_2):
    if library_1 == library_2:
        return 1.0
    else:
        tmp = sorted([library_1, library_2])
        return df.jaccard[(df.library_1==tmp[0]) & (df.library_2==tmp[1])].values[0]

mat = pd.DataFrame([[library_1, library_2] for library_1 in LIBRARIES for library_2 in LIBRARIES], columns=['library_1', 'library_2'])
mat['jaccard'] = mat.apply(lambda x: pair2jaccard(*x), axis='columns')
mat.library_1 = mat.library_1.map(lambda x: sample_info.loc[x,'label'])
mat.library_2 = mat.library_2.map(lambda x: sample_info.loc[x,'label'])
mat = mat.pivot(columns='library_2', values='jaccard', index='library_1')

def label_to_tn5(label):
    return re.match('(.*) Tn5, rep.*', label).group(1)
row_colors = None
tn5_colors = ['#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000']
tn5_colors = {tn5: color for tn5, color in zip(['0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'], tn5_colors)}
row_colors = [tn5_colors[i] for i in mat.index.map(label_to_tn5)]
p = sns.clustermap(mat, cmap='viridis', annot=False, cbar_kws={'label': 'Pairwise peak\nJaccard index'}, row_colors=row_colors)
p.savefig('{PREFIX}heatmap.pdf'.format(**locals()))
p.fig.clear()


# look at the average pairwise jaccard between tn5 concentrations...
#pd.melt(mat.reset_in)
tn5_type = pd.api.types.CategoricalDtype(categories=['0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'], ordered=True)
avg_pairwise = pd.melt(mat.reset_index(), id_vars=['library_1'])
avg_pairwise = avg_pairwise[avg_pairwise.library_1 != avg_pairwise.library_2]
avg_pairwise['tn5_1'] = avg_pairwise.library_1.map(label_to_tn5)
avg_pairwise['tn5_2'] = avg_pairwise.library_2.map(label_to_tn5)
avg_pairwise.tn5_1 = avg_pairwise.tn5_1.astype(tn5_type)
avg_pairwise.tn5_2 = avg_pairwise.tn5_2.astype(tn5_type)
avg_pairwise = avg_pairwise.groupby(by=['tn5_1', 'tn5_2'])[['value']].mean().reset_index()


# In[31]:

tmp = avg_pairwise[avg_pairwise.tn5_1==avg_pairwise.tn5_2]
tmp['tn5'] = tmp.tn5_1.map(lambda x: float(str(x).replace('X', ''))).astype(float)
plt.scatter(x=np.log2(tmp.tn5),y=tmp.value)
plt.xlabel('log2(Relative Tn5)')
plt.ylabel('Mean peak pairwise Jaccard index\nbetween experiments\nwith given Tn5 concentration')
plt.savefig('{PREFIX}tn5-vs-mean-pairwise-jaccard.pdf'.format(**locals()))
plt.clf()


# In[32]:

tmp = avg_pairwise.pivot(columns='tn5_2', values='value', index='tn5_1')
tmp.columns = tmp.columns.astype(tn5_type)
p = sns.heatmap(tmp)
p.figure.savefig('{PREFIX}mean-pairwise-jaccard-heatmap.pdf'.format(**locals()))
p.figure.clear()
