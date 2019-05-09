#!/usr/bin/env python2

import sys
import json
import gzip
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--files', nargs = '*', help = 'Names of ataqv metric files from which to extract metrics. These files will be within the data/ directory of an ataqv app.')
args = parser.parse_args()


for ataqv_data_file in args.files:

    with gzip.open(ataqv_data_file, 'rb') as f:
        line = f.readline()
        x = json.loads(line)
        tss_coverage = x[0]['tss_coverage']

        for bp in tss_coverage:
            position = int(bp[0]) - 1001
            coverage = bp[1]
            print('{}\t{}\t{}'.format(ataqv_data_file, position, coverage))
