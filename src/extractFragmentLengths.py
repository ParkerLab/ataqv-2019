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
        fragment_length_counts = x[0]['fragment_length_counts']

        fragment_lengths = [str(x) for x in sorted([int(x) for x in fragment_length_counts.keys()])]

        for fragment_length in fragment_lengths:
            count = fragment_length_counts[fragment_length][0]
            proportion = fragment_length_counts[fragment_length][1]
            print('{}\t{}\t{}\t{}'.format(ataqv_data_file, fragment_length, count, proportion))
