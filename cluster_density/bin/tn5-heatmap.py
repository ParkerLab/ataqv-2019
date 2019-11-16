#!/usr/bin/env python
# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
import sys
import os
import pybedtools as bt
import argparse
import Ataqv
import glob
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
from pygenometracks import tracks
import gffutils
import re
import itertools

parser = argparse.ArgumentParser(description='Make the public data heatmap.', add_help = True)
parser.add_argument('--sample-info', dest='sample_info', type = str, default = '', help = 'Path to public data sample info file.')
parser.add_argument('--gtf-db', dest='gtf_db', type = str, default = '', help = 'Path to the database created by gffutils.')
parser.add_argument('--heatmap-files', dest='heatmap_files', type = str, nargs = '+', help = 'List of heatmap files (*coveraged.bed).')
parser.add_argument('--ataqv-files', dest='ataqv_files', type = str, nargs = '+', help = 'List of ataqv files.')
parser.add_argument('--chrom', dest='chrom', type = str, help = 'Chromosome.')
parser.add_argument('--out', dest='out', type = str, help = 'Out file name.')
parser.add_argument('--start', dest='start', type = int, help = 'Start.')
parser.add_argument('--end', dest='end', type = int, help = 'End.')
parser.add_argument('--height', type = float, default = 5, help = 'Fig height (inches).')
parser.add_argument('--width', type = float, default = 7, help = 'Fig width (inches).')
parser.add_argument('--heatmap-and-gene-model-only', dest = 'heatmap_and_gene_model_only', action = 'store_true', default = False, help = 'Plot the heatmap and gene model (no tss enrichment or project panels)')
args = parser.parse_args()

#%matplotlib inline

SAMPLE_INFO = args.sample_info
HEATMAP_AND_GENE_MODEL_ONLY = args.heatmap_and_gene_model_only
GTF_DB = args.gtf_db
HEATMAP_FILES = args.heatmap_files
ATAQV_FILES = args.ataqv_files
CHROM = args.chrom
START = args.start
END = args.end
FIG_WIDTH = args.width
FIG_HEIGHT = args.height
OUT = args.out


# In[14]:

db = gffutils.FeatureDB(GTF_DB, keep_order=True)

# get items in the GAPDH area....
in_region = list(db.region(region=(CHROM, START, END), completely_within=False))
transcripts = [x for x in in_region if x.featuretype == 'transcript' and x.attributes.get('transcript_type')[0] == 'protein_coding']

with open('test.bed', 'w') as f:
    bed12 = '\n'.join([db.bed12(x.id, name_field='gene_name') for x in transcripts])
    f.write(bed12)


# In[16]:

# for each experiment, need the project, tss enrichment and hqaa (ataqv), and heatmap info
tn5_type = pd.api.types.CategoricalDtype(categories=['0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'])
sample_info = pd.read_csv(SAMPLE_INFO, delimiter='\t')[['library', 'tn5']].drop_duplicates().set_index('library')
sample_info.tn5 = sample_info.tn5.astype(tn5_type)
sample_info.index = sample_info.index.astype(str)

ataqv_files = []
for f in ATAQV_FILES:
    for i in sample_info.index.values:
        if i in f:
            ataqv_files.append(f)

ataqv = list(itertools.chain(*[Ataqv.load_ataqv_json(f) for f in ataqv_files]))


# In[17]:

ataqv_df = pd.DataFrame({metric: [getattr(a, metric) for a in ataqv] for metric in ['tss_enrichment', 'library', 'hqaa']}).set_index('library')
ataqv_df.index = ataqv_df.index.astype(str)
ataqv_df = ataqv_df.join(sample_info).reset_index().sort_values(['tn5', 'library']).set_index('library')
#ataqv_df = ataqv_df[ataqv_df.hqaa>=5e6]


# In[41]:

def hex2rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))


# In[63]:

tss_enrichment = ataqv_df.loc[:,['tss_enrichment']]
tn5_colors = ['#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000']
tn5_colors = {tn5: color for tn5, color in zip(['0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'], tn5_colors)}
tn5 = [hex2rgb(i) for i in ataqv_df.tn5.map(tn5_colors.get)]
tn5 = np.expand_dims(np.array(tn5), 1)


# In[24]:

# read in, reformat, and sort the heatmap...

heatmap = pd.concat([pd.read_csv(f, delimiter='\t', header=None).assign(library=re.sub(pattern=r'\.(.*)_coverage.bed', repl=r'', string=os.path.basename(f))) for f in HEATMAP_FILES])
heatmap.columns = ['chrom', 'start', 'end', 'score', 'library']
heatmap = heatmap[heatmap.library.isin(list(ataqv_df.index))]
heatmap['million_reads'] = heatmap.library.map(lambda x: ataqv_df.loc[x,'hqaa'] / 1e6)
heatmap.score = heatmap.score / heatmap.million_reads
heatmap = heatmap.drop(columns=['million_reads'])
heatmap['center'] = heatmap[['start', 'end']].apply(np.mean, axis='columns')
heatmap = heatmap.drop(columns=['chrom', 'start', 'end']).set_index('library')
heatmap = heatmap.pivot(columns='center', values='score')
heatmap = heatmap.loc[list(ataqv_df.index)]
assert(all(heatmap.index == ataqv_df.index))

# trim to the area that we will show
assert(heatmap.columns.min() <= START)
assert(heatmap.columns.max() >= END)
heatmap = heatmap.loc[:,heatmap.columns.map(lambda x: x>=START)]
heatmap = heatmap.loc[:,heatmap.columns.map(lambda x: x<=END)]

#print(heatmap.columns.min())
#print(START)
#print(heatmap.columns.max())
#print(END)


# In[25]:

heatmap = heatmap.loc[list(ataqv_df.index)]
assert(all(heatmap.index == ataqv_df.index))


# In[ ]:

N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(256/256, 1, N)
vals[:, 1] = np.linspace(0/256, 1, N)
vals[:, 2] = np.linspace(0/256, 1, N)
newcmp = mcolors.ListedColormap(np.flip(vals))


# In[61]:

# make the figure layout...
if HEATMAP_AND_GENE_MODEL_ONLY:
    fig2 = plt.figure(constrained_layout=False, figsize=(FIG_WIDTH, FIG_HEIGHT))
    spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2, top=1, bottom=0.98)
    spec2 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig2, height_ratios=[1,0.5], hspace=0, wspace=0.05, top=0.8, bottom=0)
    f2_ax3 = fig2.add_subplot(spec2[0, 0])
    f2_ax4 = fig2.add_subplot(spec2[1, 0])
    cuts_colorbar_ax = fig2.add_subplot(spec1[0, 0])

    # add gene model
    properties_dict = {'file': 'test.bed', 'color':'grey', 'labels': 'on', 'style': 'UCSC', 'fontsize': 9}
    bed = tracks.BedTrack(properties_dict)

    # axis 1 is the project
    # axis 2 is the tss enrichment
    # axis 3 is the heatmap
    # axis 4 is the gene model
    cut_norm = mcolors.Normalize(vmin=0, vmax=heatmap.max().max())
    
    f2_ax3.imshow(heatmap, aspect='auto', cmap=newcmp, extent=(heatmap.columns.min(), heatmap.columns.max(), 0, 1), norm=cut_norm)
    f2_ax3.set_xlim(START, END)
    f2_ax3.xaxis.set_ticklabels('')
    f2_ax3.yaxis.set_ticklabels('')
    f2_ax3.xaxis.set_tick_params(length=0)
    f2_ax3.yaxis.set_tick_params(length=0)

    f2_ax4.set_xlim(START, END)
    f2_ax4.set_xlabel('Chromosomal coordinate')
    f2_ax4.yaxis.set_ticklabels('')
    f2_ax4.yaxis.set_tick_params(length=0)
    bed.plot(ax=f2_ax4, chrom_region=CHROM, start_region=START, end_region=END)
    f2_ax4.ticklabel_format(axis='x', style='plain', useOffset=False)
    for t in f2_ax4.axes.get_xmajorticklabels():
        plt.setp(t, rotation=45)

    # add legends
    fig2.colorbar(mappable=cm.ScalarMappable(cmap=newcmp, norm=cut_norm), orientation='horizontal', aspect=100, cax=cuts_colorbar_ax)
    cuts_colorbar_ax.set_title('Cuts per million mapped reads')
    fig2.savefig(OUT, bbox_inches='tight', pad_inches=0)
else:
    fig2 = plt.figure(constrained_layout=False, figsize=(FIG_WIDTH, FIG_HEIGHT))
    spec1 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig2, top=1, bottom=0.98)
    spec2 = gridspec.GridSpec(ncols=3, nrows=2, figure=fig2, width_ratios=[1,1,10], height_ratios=[1,0.5], hspace=0, wspace=0.05, top=0.8, bottom=0)
    f2_ax1 = fig2.add_subplot(spec2[0, 0])
    f2_ax2 = fig2.add_subplot(spec2[0, 1])
    f2_ax3 = fig2.add_subplot(spec2[0, 2])
    f2_ax4 = fig2.add_subplot(spec2[1, 2])
    tss_enrichment_colorbar_ax = fig2.add_subplot(spec1[0, 0])
    cuts_colorbar_ax = fig2.add_subplot(spec1[0, 1])

    # add gene model
    properties_dict = {'file': 'test.bed', 'color':'grey', 'labels': 'on', 'style': 'UCSC', 'fontsize': 9}
    bed = tracks.BedTrack(properties_dict)

    # axis 1 is the project
    # axis 2 is the tss enrichment
    # axis 3 is the heatmap
    # axis 4 is the gene model

    tss_enrichment_norm = mcolors.Normalize(vmin=0, vmax=tss_enrichment.tss_enrichment.max())
    cut_norm = mcolors.Normalize(vmin=0, vmax=heatmap.max().max())

    f2_ax1.set_title('[Tn5]', fontdict={'size': 8})
    f2_ax1.imshow(tn5, aspect='auto')
    f2_ax1.xaxis.set_ticklabels('')
    f2_ax1.yaxis.set_ticklabels('')
    f2_ax1.xaxis.set_tick_params(length=0)
    f2_ax1.yaxis.set_tick_params(length=0)
    f2_ax1.set_ylabel('42 libraries')

    f2_ax2.imshow(tss_enrichment, cmap='Greys', aspect='auto', norm=tss_enrichment_norm)
    f2_ax2.set_title('TSS\nenrich.', fontdict={'size': 8})
    f2_ax2.xaxis.set_ticklabels('')
    f2_ax2.yaxis.set_ticklabels('')
    f2_ax2.xaxis.set_tick_params(length=0)
    f2_ax2.yaxis.set_tick_params(length=0)

    f2_ax3.imshow(heatmap, aspect='auto', cmap=newcmp, extent=(heatmap.columns.min(), heatmap.columns.max(), 0, 1), norm=cut_norm)
    f2_ax3.set_xlim(START, END)
    f2_ax3.xaxis.set_ticklabels('')
    f2_ax3.yaxis.set_ticklabels('')
    f2_ax3.xaxis.set_tick_params(length=0)
    f2_ax3.yaxis.set_tick_params(length=0)

    f2_ax4.set_xlim(START, END)
    f2_ax4.set_xlabel('Chromosomal coordinate')
    f2_ax4.yaxis.set_ticklabels('')
    f2_ax4.yaxis.set_tick_params(length=0)
    bed.plot(ax=f2_ax4, chrom_region=CHROM, start_region=START, end_region=END)
    f2_ax4.ticklabel_format(axis='x', style='plain', useOffset=False)
    for t in f2_ax4.axes.get_xmajorticklabels():
        plt.setp(t, rotation=45)


    # add legends
    fig2.colorbar(mappable=cm.ScalarMappable(cmap='Greys', norm=tss_enrichment_norm), orientation='horizontal', aspect=100, cax=tss_enrichment_colorbar_ax)
    fig2.colorbar(mappable=cm.ScalarMappable(cmap=newcmp, norm=cut_norm), orientation='horizontal', aspect=100, cax=cuts_colorbar_ax)
    tss_enrichment_colorbar_ax.set_title('TSS enrichment')
    cuts_colorbar_ax.set_title('Cuts per million mapped reads')
    fig2.savefig(OUT, bbox_inches='tight', pad_inches=0)


# In[ ]:



