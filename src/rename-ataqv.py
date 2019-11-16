#!/usr/bin/env python 

import sys
import json
import gzip
import argparse

parser = argparse.ArgumentParser(description='Update the name in an ataqv metrics file.', add_help = True)
parser.add_argument('metrics_file', type = str, help = 'ataqv metric file.')
parser.add_argument('new_name', type = str,  help = 'New name.')
args = parser.parse_args()

with gzip.open(args.metrics_file, 'rb') as f:
        line = f.read()
        d = json.loads(line)[0]
        d['metrics']['library']['sample'] = args.new_name
        print(json.dumps([d], indent = 4, sort_keys = True))
