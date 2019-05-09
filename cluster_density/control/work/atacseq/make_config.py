#!/usr/bin/env python

import os
import sys
import json
import csv
import re

ROOT = sys.argv[1]
RESULTS = sys.argv[2]
FASTQ_DIR = sys.argv[3]
FASTQ_FILES = sys.argv[4:]
GENOME = 'hg19'
DATA = os.path.join(ROOT, 'data')


def fastq_to_readgroup(fastq):
    return re.match('^(.*)\.[12]\.fastq\.gz$', os.path.basename(fastq)).group(1)


def readgroup_to_library(readgroup):
    return re.match('^(\d+)___L\d+___\d+$', readgroup).group(1)


def readgroup_to_sequencing_run(readgroup):
    return re.match('^\d+___L\d+___(\d+)$', readgroup).group(1)


def get_fastq_files(readgroup):
    tmp = os.path.join(FASTQ_DIR, readgroup)
    return [tmp + '.{}.fastq.gz'.format(x) for x in [1, 2]]


libraries = dict()

for fastq in FASTQ_FILES:
    readgroup = fastq_to_readgroup(fastq)
    library = readgroup_to_library(readgroup)
    sequencing_run = readgroup_to_sequencing_run(readgroup)
    key = '{library}___{sequencing_run}'.format(**locals())

    if key not in libraries:
        libraries[key] = {
            'genome': GENOME,
            'readgroups': dict()
        }

    libraries[key]['readgroups'][readgroup] = get_fastq_files(readgroup)


MAPPABILITY_DIR = os.path.join(DATA, 'mappability')

CONFIG = {
    'libraries': libraries,
    'results': RESULTS,
    'bwa_index': {
        GENOME: os.path.join(DATA, 'bwa', GENOME, GENOME)
    },
    'whitelist': {
        GENOME: [os.path.join(MAPPABILITY_DIR, x) for x in os.listdir(MAPPABILITY_DIR) if re.match('^{}.whitelist'.format(GENOME), x)]
    },
    'blacklist': {
        GENOME: [os.path.join(MAPPABILITY_DIR, x) for x in os.listdir(MAPPABILITY_DIR) if re.match('^{}.blacklist'.format(GENOME), x)]
    },
    'tss': {
        GENOME: os.path.join(ROOT, 'src', 'ataqv', 'data', 'tss', '{}.tss.refseq.bed.gz'.format(GENOME))
    }
}

print(json.dumps(CONFIG, sort_keys = True, indent = 4))
