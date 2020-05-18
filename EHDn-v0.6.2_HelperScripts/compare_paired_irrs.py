#!/usr/bin/env python3

#
# Expansion Hunter Denovo
# Copyright (c) 2017 Illumina, Inc.
#
# Author: Egor Dolzhenko <edolzhenko@illumina.com>,
#         Sai Chen <schen6@illumina.com>
# Concept: Michael Eberle <meberle@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import argparse
import logging
import sys

from core import common


def load_parameters():
    '''Capture command-line parameters.

    '''
    parser = argparse.ArgumentParser(
        description='compare counts of in-repeat read pairs')
    parser.add_argument(
        '--manifest',
        help='TSV file with id, case/control status, and path of each sample',
        required=True)
    parser.add_argument(
        '--combinedCounts',
        help='file with combined counts of anchored in-repeat reads',
        required=True)
    parser.add_argument(
        '--pirRegions',
        help='BED file with results',
        required=True)
    parser.add_argument(
        '--minCount',
        help='minimum number reads in a region for downstream analysis',
        default=5,
        type=int)

    args = parser.parse_args()

    parameter_encoding = ' '.join(sys.argv[1:])
    logging.info('Starting with these parameters: %s', parameter_encoding)

    return {'manifest_path': args.manifest, 'counts_path': args.combinedCounts,
            'min_count': args.minCount, 'output_path': args.pirRegions}


def generate_count_table(combined_counts):
    count_table = []
    for unit, rec in combined_counts.items():
        if 'IrrPairCounts' not in rec:
            continue
        sample_counts = rec['IrrPairCounts']
        table_row = {'unit': unit, 'sample_counts': sample_counts}
        count_table.append(table_row)

    logging.info('Loaded %i repeat units', len(count_table))
    return count_table


def output_results(count_table, output_path):
    with open(output_path, 'w') as output_file:
        for row in count_table:
            unit = row['unit']
            pvalue, bonf_pvalue = row['pvalue'], row['bonf_pvalue']

            sample_counts = row['sample_counts']
            encoded_counts = ['{}:{}'.format(s, c)
                              for s, c in sample_counts.items()]
            encoded_counts = ','.join(encoded_counts)
            print(unit, pvalue, bonf_pvalue, encoded_counts,
                  sep='\t', file=output_file)


def main():
    common.init_logger()
    parameters = load_parameters()
    combined_counts = common.load_combined_json(parameters['counts_path'])
    samples = common.load_manifest(parameters['manifest_path'])
    count_table = generate_count_table(combined_counts)
    count_table = common.filter_counts(parameters['min_count'], count_table)
    sample_status = common.extract_case_control_assignments(samples)
    common.compare_counts(sample_status, count_table)
    common.correct_pvalues(count_table)
    output_results(count_table, parameters['output_path'])
    logging.info('Done')


if __name__ == '__main__':
    main()
