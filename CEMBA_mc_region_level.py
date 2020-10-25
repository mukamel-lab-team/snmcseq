#!/usr/bin/python3
"""
# Read in allc tables (with tabix open) and generate summarized mc levels within gene bodies
"""

from __init__ import *
# import pandas as pd
# import numpy as np
# import os
# import time
import tabix
import subprocess as sp
import argparse

import snmcseq_utils
from snmcseq_utils import create_logger

def mc_region_level_worker(allc_file, output_file, bed_file, chr_prefix=False, bed_file_name_column=False,
    contexts=CONTEXTS, compress=True, cap=2, species=SPECIES):
    """
    allc_file 
    bed file:
    """
    # logger = create_logger()

    chromosomes = snmcseq_utils.get_chromosomes(species)  

    if bed_file_name_column:
        columns = ['chr', 'start', 'end', 'name']
        df_gtf = pd.read_table(bed_file, header=None, 
                names=columns, usecols=[0, 1, 2, 3], dtype={'chr': object, 'name': object})
    else:
        columns = ['chr', 'start', 'end']
        df_gtf = pd.read_table(bed_file, header=None, 
                names=columns, usecols=[0, 1, 2], dtype={'chr': object})



    df_gtf.chr = df_gtf.chr.apply(lambda x: x[len('chr'):] if x.startswith('chr') else x)
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    logging.info("mc_region_level processing: {} {}".format(allc_file, contexts))


    outfile = open(output_file, "w")
    for context in contexts:
        columns += ['m{}'.format(context), context]

    outfile.write('\t'.join(columns)+'\n')

    allc = tabix.open(allc_file)
    for i, row in df_gtf.iterrows():
        if bed_file_name_column:
            row_out = [str(row.chr), str(row.start), str(row.end), str(row['name'])]
        else:
            row_out = [str(row.chr), str(row.start), str(row.end)]

        if chr_prefix: 
            records = list(allc.query('chr'+str(row['chr']), row['start'], row['end']))
        else:
            records = list(allc.query(row['chr'], row['start'], row['end']))

        for context in contexts:
            mc, c = snmcseq_utils.tabix_summary(records, context=context, cap=cap) # remove sites with total_c > 2
            row_out += [str(mc), str(c)] 

        outfile.write('\t'.join(row_out)+'\n')

    outfile.close()

    logging.info("Done with mc_region_level processing: {} {}\n Saving results to {}".format(allc_file, contexts, output_file))

    if compress:
        snmcseq_utils.compress(output_file, suffix='gz') # Fangming 12/3/2018 

    return 0

def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc", help="allc file path", required=True)
    parser.add_argument("-o", "--output_file", help="output file path", required=True)
    parser.add_argument("-b", "--bed_file", help="bed file", required=True)
    parser.add_argument("-c", "--contexts", help="list of contexts: CH/CG/...", nargs='+', default=CONTEXTS)
    # parser.add_argument("-f", "--overwrite", 
    #     action='store_true',
    #     help="overwrite a file if it exists")
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()

    log = create_logger()

    # allc_file = '/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3C_171206/allc/allc_171213_CEMBA_mm_P56_P63_3C_MOp_CEMBA171206_3C_1_CEMBA171206_3C_3_H4_AD002_indexed.tsv.bgz'
    # output_file = 'test.tsv'
    # bed_file = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/CEMBA_3C_171206/CEMBA_3C_171206_merged.bed'

    # mc_region_level_worker(allc_file, output_file, bed_file, contexts=CONTEXTS)

    mc_region_level_worker(args.input_allc, args.output_file, args.bed_file, contexts=['CH', 'CG'])


