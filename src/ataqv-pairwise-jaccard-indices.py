
# coding: utf-8

# In[76]:

import pybedtools as bt
import numpy as np
import pandas as pd
import sys
import os
from matplotlib import pyplot as plt
import seaborn as sns
import itertools
import re

# %matplotlib inline

SAMPLE_INFO = sys.argv[1]
OUT = sys.argv[2]
PEAK_FILES = sorted(sys.argv[3:])


sample_info = pd.read_csv(SAMPLE_INFO, delimiter='\t')[['tn5', 'replicate', 'library']]
sample_info.library = sample_info.library.map(str)
sample_info = sample_info.set_index('library')
sample_info['label'] = sample_info[['tn5', 'replicate']].apply(lambda x: '{} Tn5, {}'.format(*x), axis='columns')


# In[57]:

def peak_file_to_library(f):
    return re.match('^(.*)_peaks.broadPeak.*', os.path.basename(f)).group(1)


# In[59]:

def jaccard(peaks_1, peaks_2):
    # should be bedtools objects
    # merge to create a set of master peaks
    merge = peaks_1.cat(peaks_2).merge()
    in_peaks_1 = merge.intersect(peaks_1, u=True)
    in_peaks_2 = merge.intersect(peaks_2, u=True)
    number_shared = in_peaks_1.intersect(in_peaks_2, u=True).count()
    #number_specific_to_1 = in_peaks_1.intersect(in_peaks_2, v=True).count()
    #number_specific_to_2 = in_peaks_2.intersect(in_peaks_1, v=True).count()
    jaccard_index = number_shared / merge.count()
    return jaccard_index


# In[60]:

# calculate jaccard index for each pair...
peak_file_pairs = pd.DataFrame([[x, y] for x, y in list(itertools.combinations(PEAK_FILES, 2))], columns=['file_1', 'file_2'])
peak_file_pairs['jaccard'] = peak_file_pairs.apply(lambda x: jaccard(*[bt.BedTool(i) for i in x]), axis='columns')


# In[61]:

peak_file_pairs['library_1'] = peak_file_pairs.file_1.map(peak_file_to_library)
peak_file_pairs['library_2'] = peak_file_pairs.file_2.map(peak_file_to_library)
peak_file_pairs = peak_file_pairs.drop(columns=['file_1', 'file_2'])


# In[66]:

def pair2jaccard(library_1, library_2):
    if library_1 == library_2:
        return 1.0
    else:
        tmp = sorted([library_1, library_2])
        return peak_file_pairs.jaccard[(peak_file_pairs.library_1==tmp[0]) & (peak_file_pairs.library_2==tmp[1])].values[0]


# In[84]:

mat = pd.DataFrame([[peak_file_to_library(x), peak_file_to_library(y)] for x in PEAK_FILES for y in PEAK_FILES], columns=['library_1', 'library_2'])
mat['jaccard'] = mat.apply(lambda x: pair2jaccard(*x), axis='columns')
mat.library_1 = mat.library_1.map(lambda x: sample_info.loc[x,'label'])
mat.library_2 = mat.library_2.map(lambda x: sample_info.loc[x,'label'])
mat = mat.pivot(columns='library_2', values='jaccard', index='library_1')


# In[89]:

p = sns.clustermap(mat, cmap='viridis', annot=True, cbar_kws={'label': 'Pairwise peak\nJaccard index'})
p.savefig(OUT)


# In[ ]:



