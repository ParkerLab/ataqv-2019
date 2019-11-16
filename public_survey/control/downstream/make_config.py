import os
import sys
import json
import csv
import glob

ROOT, ATAC_DIR, sample_info, RESULTS = sys.argv[1:] 

config = dict()
config['src'] = os.path.join(ROOT, 'public_survey', 'bin')
config['sample_info'] = sample_info
config['libraries'] = dict()
config['results'] = RESULTS
config['gtf'] = {'hg19': os.path.join(ROOT, 'data', 'gtf', 'hg19.gtf.gz')}
config['chrom_sizes'] = {genome: os.path.join(ROOT, 'data', 'chrom_sizes', genome + '.chrom.sizes') for genome in ['hg19', 'mm9']}
config['tss'] = {genome: os.path.join(ROOT, 'data', 'tss', genome + '.tss.bed.gz') for genome in ['hg19', 'mm9']}
config['blacklist'] = {genome: glob.glob(os.path.join(ROOT, 'data', 'mappability', genome + '.*')) for genome in ['hg19', 'mm9']}

config['atacseq'] = {
                #'ataqv': os.path.join(ATAC_DIR, 'ataqv'),
                #'ataqv_readgroups': os.path.join(ATAC_DIR, 'ataqv.readgroups'),
                'peaks': os.path.join(ATAC_DIR, 'macs2'),
                'md': os.path.join(ATAC_DIR, 'mark_duplicates'),
                'pruned': os.path.join(ATAC_DIR, 'prune')
            }

with open(sample_info, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        genome = 'hg19' if line['organism'] == 'Homo sapiens' else 'mm9'
        x = 'single-cell' if line['is_single_cell'] == 'TRUE' else 'bulk'
        if line['experiment'] not in config['libraries']:
            config['libraries'][line['experiment']] = {'genome': genome, 'project': line['project'], 'class': x, 'runs': []}
        config['libraries'][line['experiment']]['runs'].append(line['run'])

print(json.dumps(config, indent = 4, sort_keys = True))
