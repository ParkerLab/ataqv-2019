#!/usr/bin/env python
# coding: utf-8

# In[71]:

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import sys
import scipy.stats as stats
import statsmodels.api as sm
import seaborn as sns

RESULTS_WITH_NO_COVARIATE, RESULTS_WITH_COVARIATE = sys.argv[1:]

def expected_p_values(observed_p_values):
    return list(stats.rankdata(observed_p_values, method='ordinal')/ (len(observed_p_values)+1))


no_covariate = pd.read_csv(RESULTS_WITH_NO_COVARIATE, delimiter='\t').set_index('peak')
no_covariate['covariates'] = 'none'
no_covariate['expected'] = expected_p_values(no_covariate.pvalue)
no_covariate['x'] = -1*np.log10(no_covariate.expected)
no_covariate['y'] = -1*np.log10(no_covariate.pvalue)

mfl_covariate = pd.read_csv(RESULTS_WITH_COVARIATE, delimiter='\t').set_index('peak')
mfl_covariate['covariates'] = 'median fragment length'
mfl_covariate['expected'] = expected_p_values(mfl_covariate.pvalue)
mfl_covariate['x'] = -1*np.log10(mfl_covariate.expected)
mfl_covariate['y'] = -1*np.log10(mfl_covariate.pvalue)

both = pd.concat([no_covariate, mfl_covariate])


no_covariate.pvalue.hist()
plt.savefig('no-covariate-p-value-distribution.pdf')
plt.clf()

mfl_covariate.pvalue.hist()
plt.savefig('median-fragment-length-covariate-p-value-distribution.pdf')
plt.clf()

both['significant'] = both.padj.map(lambda x: x <= 0.05)
both['direction'] = np.sign(both.log2FoldChange)
# both.groupby(by=['significant', 'covariates', 'direction']).count()

p = sns.relplot(x='x', y='y', hue='covariates', alpha=0.3, data=both)
plt.sca(p.ax)
for_diagonal_line = min([both.x.max(), both.y.max()])
plt.plot([0, for_diagonal_line], [0, for_diagonal_line], c='black')
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.savefig('su-et-al-qqplot.pdf')
