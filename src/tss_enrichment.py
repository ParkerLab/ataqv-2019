#!/usr/bin/env python

import sys
import os
import pybedtools
from pybedtools.featurefuncs import five_prime, three_prime
import argparse
import random
import os

parser = argparse.ArgumentParser(description='Calculate TSS enrichment from a bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file, sorted by read name.')
parser.add_argument('tss', type = str, help = 'BED file of TSS.')
parser.add_argument('genome', type = str, help = 'Chromosome size file (bedtools genome file)')
parser.add_argument('--method', type = str, default = 'ENCODE', help = 'Method to use (ENCODE or ataqv)')
parser.add_argument('--read-length', dest = 'read_length', type = int, default = 100, help = 'Read length (used only for ENCODE method). Default = 100')

args = parser.parse_args()

tmp_rand = random.randint(1, 1000000)

def get_fragments(bam):
    fragments = pybedtools.BedTool(bam).bam_to_bed(bedpe = True).cut([0, 1, 5]).sort()
    return(fragments)


def get_cutsite_centered_reads(bam, read_length, genome):
    # require that each shifted read lie completely within the chromosome coordinates
    # Therefore, filter out cutsites that would break this.
    
    chromosome_ends = ''
    with open(genome, 'r') as g:
        for line in g:
            chrom, length = line.rstrip().split('\t')
            length = int(length)
            chromosome_ends += '{}\t{}\t{}\n'.format(chrom, 0, int(0+read_length/2))
            chromosome_ends += '{}\t{}\t{}\n'.format(chrom, int(length-read_length/2), length)

    # names of tmp files
    chromosome_ends_file = 'chromosome_ends.{}.bed'.format(tmp_rand)
    fragments_file = 'fragments.{}.bed'.format(tmp_rand)
    five_prime_file = 'five_prime.{}.bed'.format(tmp_rand)
    three_prime_file = 'three_prime.{}.bed'.format(tmp_rand)
    
    pybedtools.BedTool(chromosome_ends, from_string = True).saveas(chromosome_ends_file)

    pybedtools.BedTool(bam).bam_to_bed(bedpe = True).cut([0, 1, 5]).saveas(fragments_file)
    pybedtools.BedTool(fragments_file).each(five_prime, upstream = 0, downstream = 1).saveas(five_prime_file)
    pybedtools.BedTool(fragments_file).each(three_prime, upstream = 0, downstream = 1).saveas(three_prime_file)
    cutsite_1 = pybedtools.BedTool(five_prime_file).intersect(b = chromosome_ends_file, v = True).slop(b = int(read_length / 2), g = genome)
    cutsite_2 = pybedtools.BedTool(three_prime_file).intersect(b = chromosome_ends_file, v = True).slop(b = int(read_length / 2), g = genome)
    fragments = cutsite_1.cat(cutsite_2, postmerge = False).sort()
    os.remove(chromosome_ends_file)
    os.remove(fragments_file)
    os.remove(five_prime_file)
    os.remove(three_prime_file)

    return fragments


sys.stderr.write('Fetching reads...\n')

reads = None

if args.method == 'ataqv':
    reads = get_fragments(args.bam)
elif args.method == 'ENCODE':
    reads = get_cutsite_centered_reads(args.bam, args.read_length, args.genome)
else:
    sys.stderr.write('Unrecognized method: {}'.format(args.method))
    assert(False)

sys.stderr.write('Reading in TSS...\n')
sys.stderr.flush()

tss = pybedtools.BedTool(args.tss).slop(b = 1000, g = args.genome).sort()

sys.stderr.write('Getting coverages...\n')
sys.stderr.flush()

c = tss.coverage(reads, d = True, sorted = True)

sys.stderr.write('Calculating enrichments...\n')
sys.stderr.flush()

sums = {}
for i in range(-1000, 1001):
    sums[i] = 0

for i in c:
    strand = i[5]
    signal = int(i[7])

    position = int(i[6]) - 1001
    if strand == '-':
        position = 1001 - int(i[6])
    sums[position] += signal

baseline = 0

for i in range(-1000, -900):
    baseline += sums[i]

for i in range(901, 1001):
    baseline += sums[i]

baseline = float(baseline) / 200

for i in range(-1000, 1001):
    print('{}\t{}'.format(i, sums[i] / baseline))
