#!/usr/bin/env python

# subsample a set of bam files to the same depth
# assumes that any filtering etc has already been done -- i.e.,
# the bams should be subsetted to the minimum value of ...

import pysam
import sys
import os
import argparse
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Subsample a set of bam files.')
parser.add_argument('--suffix', type=str, default = 'subsampled', help='How to name the output files. XXXXX.bam -> XXXXX.suffix.bam. Default: subsampled')
parser.add_argument('bam', nargs='+', type=str, help='The bams to subsample.')
parser.add_argument('--number-reads', dest = 'number_reads', type=int, nargs='?', help='The (approximate) number of reads to subsample to. Either this argument or --same-depth must be given.')
parser.add_argument('--same-depth', dest = 'same_depth', action='store_true', help='If this flag is passed, subsample the files to the same depth (= the depth of the most shallow file). Ignored if --number-reads is passed.')
parser.add_argument('--seed', type=int, default=872674, help='Seed for the sampler.')
args = parser.parse_args()

if args.number_reads and args.same_depth:
    logging.warning("--number-reads was passed; ignoring --same-depth")
    args.same_depth = False

if args.number_reads is None and not args.same_depth:
    logging.error("--number-reads or --same-depth must be passed")
    sys.exit(1)

# Determine the number of QC-passed reads in each of the bam files
logging.info("Determining the number of QC-passed reads in each file")

flagstats = [pysam.flagstat(x).split('\n')[0] for x in args.bam]
number_reads = [int(re.match('\A(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)\Z', i).group(1)) for i in flagstats]

for index, bam in enumerate(args.bam):
    logging.info("Found {} reads in {}".format(number_reads[index], bam))


# Now perform the subsampling
for index, bam in enumerate(args.bam):

    bam_prefix = re.match('\A(.*).bam\Z', os.path.basename(bam)).group(1)
    subsample_fraction = None

    if args.number_reads is not None:
        subsample_fraction = float(args.number_reads) / number_reads[index]
    else:
        subsample_fraction = min(number_reads) / float(number_reads[index])

    if args.number_reads:
        if subsample_fraction >= 1:
            logging.info("Keeping all reads from {}, since the number of reads desired is greater than or equal to the number of reads in the original bam".format(bam))
        else:
            logging.info("Subsampling {} to approximately {} reads".format(bam, args.number_reads))
    else:
        logging.info("Subsampling {} to approximately {} reads".format(bam, min(number_reads)))


    if subsample_fraction >= 1:
        pysam.view('-h', '-b', '-o', '{}.{}.bam'.format(bam_prefix, args.suffix), bam, catch_stdout = False)
    else:
        subsample_fraction = '{0:f}'.format(subsample_fraction)
        subsample_fraction = re.match('\A[10]\.(.*)\Z', subsample_fraction).group(1)

        pysam.view('-h', '-b', '-s', '{}.{}'.format(args.seed, subsample_fraction), '-o', '{}.{}.bam'.format(bam_prefix, args.suffix), bam, catch_stdout = False)
