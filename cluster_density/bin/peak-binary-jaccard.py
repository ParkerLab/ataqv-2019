#!/usr/bin/env python
# coding: utf-8

# In[2]:

import pybedtools as bt
import numpy as np
import pandas as pd
import sys
import os
import itertools
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

OUT = sys.argv[1]
PEAK_FILES = sorted(sys.argv[2:])

def jaccard(file_1, file_2):
    logging.info('Calculating binary Jaccard between {file_1} and {file_2}'.format(**locals()))
    peaks_1 = bt.BedTool(file_1)
    peaks_2 = bt.BedTool(file_2)
    merge = peaks_1.cat(peaks_2).merge()
    in_peaks_1 = merge.intersect(peaks_1, u=True)
    in_peaks_2 = merge.intersect(peaks_2, u=True)
    number_shared = in_peaks_1.intersect(in_peaks_2, u=True).count()
    jaccard_index = number_shared / merge.count()
    bt.helpers.cleanup()
    return jaccard_index

# calculate jaccard index for each pair...
peak_file_pairs = pd.DataFrame([[x, y] for x, y in list(itertools.combinations(PEAK_FILES, 2))], columns=['file_1', 'file_2'])
peak_file_pairs['jaccard'] = peak_file_pairs.apply(lambda x: jaccard(*x), axis='columns')

peak_file_pairs.to_csv(OUT, sep='\t', index=False)
