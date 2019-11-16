
# coding: utf-8

# In[39]:

import os
import sys
import glob
import pandas as pd
import numpy as np
import pybedtools as bt
import argparse
import re

parser = argparse.ArgumentParser(description='Setup the differential peak analysis calling and print the commands to run.', add_help = True)
parser.add_argument('--min-samples', dest='min_samples', type = int, default=1,  help = 'Minimum number of samples that must support a peak.')
parser.add_argument('--fdr', type = float, default=0.01,  help = 'FDR threshold.')
parser.add_argument('--peak-files', dest='peak_files', nargs='+', type = str, help = 'Peak files (blacklist-filtered).')

args = parser.parse_args()


def bedcat(beds, postmerge=False):
    if len(beds) == 1:
        return beds
    else:
        return beds[0].cat(*beds[1:], postmerge=postmerge)


def truncate_name(feature):
    feature.name = re.sub(r'(.*)_peak.*', r'\1', feature.name)
    return feature


# In[41]:

bed = bedcat([bt.BedTool(f) for f in args.peak_files])


# In[44]:

# fdr filter
bed = bed.filter(lambda x: float(x[8]) >= -1*np.log10(args.fdr)).saveas()


# In[48]:

# change the peak names...
bed = bed.each(truncate_name).saveas()


# In[50]:

# merge
potential_master_peaks = bed.sort().merge(c=4, o='count_distinct')


# In[58]:

master_peaks = potential_master_peaks.filter(lambda x: int(x.name) >= args.min_samples).cut([0,1,2]).saveas()
for p in master_peaks:
    print(str(p).rstrip())

