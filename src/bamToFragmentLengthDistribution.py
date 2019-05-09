#!/usr/bin/env python2

import sys
import json
import gzip
import numpy
import argparse

import pysam

parser = argparse.ArgumentParser()
parser.add_argument('bam', nargs = '*', help = 'List of bam files to analyze')
parser.add_argument('--mapq', default = 30, type = int, help = 'Only use reads with mapping quality at least this high (default: 30)')
args = parser.parse_args()


def find_index(l, i):
    """Input: [3, 2, 6], 4
    Output: 1
    Logic:
    [3, 2, 6] --> [0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2]
    4th value in that list is 1
    """

    assert(sum(l) >= i)

    running_sum = 0
    current_index = 0

    for j in l:
        if running_sum + j < i:
            running_sum += j
            current_index += 1
        else:
            return current_index


def median_from_dict(d):
    # d[value] --> count
    # this is useful for e.g. finding the median fragment length, since that is stored
    # as dict[fragment_length] -> count

    keys = sorted(d.keys())
    values = [d[i] for i in keys]

    assert(sum(values) > 0)

    # remove 0 values...
    while 0 in values:
        i = values.index(0)
        del keys[i]
        del values[i]

    count = sum(values)

    middle = None

    if count % 2 == 0:
        middle = [int(count / 2) - 1, int(count / 2)]  # we want these indices, when the values are rep()'ed and ranked
    else:
        middle = [int(count) / 2]

    indices = [find_index(values, i + 1) for i in middle]
    corresponding_keys = [keys[i] for i in indices]

    return numpy.mean(corresponding_keys)


def mode_from_dict(d):
    # d[value] --> count

    modes = []
    max_val = max(d.values())

    for key, val in d.items():
        if val == max_val:
            modes.append(key)

    return modes


def median(l):
    # l should be a list of numeric values
    return numpy.median(numpy.array(l))


def variance_from_dict(d, mu):
    numerator = 0

    for key, value in d.items():
        numerator += ( value * ((key - mu)**2) )

    return numerator / sum(d.values())


for bam in args.bam:

    insert_sizes = {}

    f = pysam.AlignmentFile(bam, 'rb')

    for read in f.fetch(until_eof=True):
        if read.mapping_quality >= args.mapq:
            template_length = read.template_length
            if template_length < 0:
                continue
            else:
                if template_length not in insert_sizes:
                    insert_sizes[template_length] = 0
                insert_sizes[template_length] += 1

    # now calculate the statistics
    # first, average insert size
    keys = insert_sizes.keys()
    vals = [insert_sizes[i] for i in keys]
    mu = numpy.average(numpy.array(keys), weights=numpy.array(vals))
    print('{}\tmean_insert_size\t{}'.format(bam, mu))
    print('{}\tmedian_insert_size\t{}'.format(bam, median_from_dict(insert_sizes)))
    print('{}\tinsert_size_mode\t{}'.format(bam, mode_from_dict(insert_sizes)))
    print('{}\tinsert_size_sd\t{}'.format(bam, variance_from_dict(insert_sizes, mu)**(0.5)))

    for key, val in sorted(insert_sizes.items()):
        print('{}\tinsert_size_hist\t{}:{}'.format(bam, key, val))
