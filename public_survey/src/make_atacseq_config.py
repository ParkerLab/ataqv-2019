#!/usr/bin/env python

import os
import sys
import csv
import json
import argparse

parser = argparse.ArgumentParser(description='Prepare a config file for ATAC-seq', add_help = True)
parser.add_argument('sample_info', type = str,  help = 'Path to sample information tsv file')
parser.add_argument('results', type = str, help = 'Path to the desired results directory.')
parser.add_argument('fastq_dir', type = str, help = 'Path to the fastq file directory.')
parser.add_argument('--run-level', dest = 'run_level', action = 'store_true', default = False, help = 'If flag is passed process runs (SRR* accessions) rather than libraries (SRX* accessions).')

args = parser.parse_args()

CONFIG = dict()
LIBRARIES = dict()

with open(args.sample_info, 'r') as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for line in reader:
        run_accession = line['run']
        library_accession = line['experiment']
        identifier = run_accession if args.run_level else library_accession
        if identifier not in LIBRARIES:
            LIBRARIES[identifier] = {
                    'genome': 'hg19' if line['organism'] == 'Homo sapiens' else 'mm9',
                    'readgroups' : dict(),
                    'is_single_cell': line['is_single_cell'],
                    'title': line['title']
                }
        LIBRARIES[identifier]['readgroups'][run_accession] = [os.path.join(args.fastq_dir, run_accession + '.{x}.fastq.gz'.format(**locals())) for x in [1,2]]


CONFIG['results'] = args.results
CONFIG['libraries'] = LIBRARIES
CONFIG['bwa_index'] = {
        'hg19': '/lab/work/porchard/data/bwa/hg19/hg19',
        'mm9': '/lab/work/porchard/data/bwa/mm9/mm9'
    }
CONFIG['tss'] = {
        'hg19': '/home/porchard/github/ataqv/data/tss/hg19.tss.refseq.housekeeping.ortho.bed.gz',
        'mm9': '/home/porchard/github/ataqv/data/tss/mm9.tss.refseq.housekeeping.ortho.bed.gz'
    }
CONFIG['blacklist'] = {
        'hg19': ['/lab/data/reference/human/hg19/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz', '/lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz'],
        'mm9': ['/lab/data/reference/mouse/mm9/annot/mm9-blacklist.bed.gz']
    }

print(json.dumps(CONFIG, sort_keys = True, indent = 4))
