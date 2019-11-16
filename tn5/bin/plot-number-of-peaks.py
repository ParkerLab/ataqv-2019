#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pybedtools as bt
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import seaborn as sns
import argparse
import re
import logging
import numpy as np


parser = argparse.ArgumentParser(description='Analyze peak calling results.', add_help = True)
parser.add_argument('--sample-info', dest = 'sample_info', type = str, help = 'Sample info file.')
parser.add_argument('--peaks', type = str, nargs='+', help = 'Set of .broadPeak files. Name: {library}_peaks.broadPeak*')
parser.add_argument('--out', type = str, help = 'Name of output PDF.')

args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

sample_info = pd.read_csv(args.sample_info, delimiter='\t')
sample_info.library = sample_info.library.map(str)
sample_info.index = sample_info.library
sample_info.tn5 = sample_info.tn5.map(lambda x: np.log2(float(x.replace('X', ''))))

def peak_file_to_library(f):
    FILENAME_RE = '^(.*)_peaks.broadPeak.*'
    return re.match(FILENAME_RE, os.path.basename(f)).group(1)


def peak_to_library(peak_name):
    return re.match('^(.*)_peak_\d+').group(1)


peaks = {peak_file_to_library(f): bt.BedTool(f) for f in args.peaks}


# first, just plot the number of peaks per library...
# for each library, determine the number of total peaks and the number of peaks in each chromatin state???
# filter down to just the libraries of interest...
sample_info = sample_info[sample_info.library.isin(list(peaks.keys()))]
sample_info['number_peaks'] = sample_info.library.map(lambda x: peaks[x].count())


# In[5]:

# map replicates to colors...
colors = {f'rep{replicate}': color for replicate, color in zip(range(1, 7), ['#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e'])}
for i in list(colors.keys()):
    if i not in sample_info.replicate.values:
        colors.pop(i)


# In[8]:

# plot the number of total peaks per experiment
# plot just the high cluster density experiment if examining the cluster density data...
for replicate in sample_info.replicate.unique():
    tmp = sample_info[sample_info.replicate==replicate]
    if 'cluster_density' in sample_info.columns:
        tmp = tmp[tmp.cluster_density=='high']
    plt.plot('tn5', 'number_peaks', color=colors[replicate], data=tmp)
    plt.scatter('tn5', 'number_peaks', color=colors[replicate], data=tmp)

plt.ylim([sample_info.number_peaks.min()*0.9, sample_info.number_peaks.max()*1.1])
plt.xlabel('log2(Relative [Tn5])')
plt.ylabel('Number of peaks')
plt.legend(handles=[mpatches.Patch(color=color, label=label) for color, label in [(color, label) for label, color in sorted(colors.items())]])
#plt.title('Number of total peaks')
plt.savefig(args.out)

