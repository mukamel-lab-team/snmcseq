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

def mc_region_level_worker(allc_file, output_file, bed_file,
    contexts=CONTEXTS, compress=True):
    """
    allc_file 
    bed file:
    """
    # logger = create_logger()

    chromosomes = snmcseq_utils.get_mouse_chromosomes() # specific to MOUSE 

    df_gtf = pd.read_table(bed_file, header=None, 
            names=['chr', 'start', 'end'], usecols=[0, 1, 2], dtype={'chr': object})
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[len('chr'):] if x.startswith('chr') else x)
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    logging.info("mc_region_level processing: {} {}".format(allc_file, contexts))


    outfile = open(output_file, "w")
    columns = ['chr', 'start', 'end']
    for context in contexts:
        columns += ['m{}'.format(context), context]

    outfile.write('\t'.join(columns)+'\n')

    allc = tabix.open(allc_file)
    for i, row in df_gtf.iterrows():
        row_out = [str(row.chr), str(row.start), str(row.end)]
        records = list(allc.query(row.chr, row.start, row.end))
        for context in contexts:
            mc, c = snmcseq_utils.tabix_summary(records, context=context, cap=2) # remove sites with total_c > 2
            row_out += [str(mc), str(c)] 

        outfile.write('\t'.join(row_out)+'\n')

    outfile.close()

    logging.info("Done with mc_region_level processing: {} {}\n Saving results to {}".format(allc_file, contexts, output_file))

    if compress:
        snmcseq_utils.compress(output_file) 

    return 0


# def mc_gene_level(allc_file,
#     contexts=CONTEXTS, 
#     genebody=GENEBODY,
#     convention='CEMBA', 
#     overwrite=False):

#     """
#     set up conventions for output_file
#     """

#     if convention=='CEMBA':
#         CEMBA_DATASETS = '/cndd/Public_Datasets/CEMBA/Datasets'
#         allc_file = os.path.abspath(allc_file)
#         assert allc_file[:len(CEMBA_DATASETS)] == CEMBA_DATASETS

#         dataset, *dis, allc_basename = allc_file[len(CEMBA_DATASETS)+1:].split('/')


#         sample = allc_basename[len('allc_'):-len('.tsv.bgz')] 

#         output_dir = "{}/{}/gene_level".format(CEMBA_DATASETS, dataset)
#         output_file = "{}/genebody_{}.tsv".format(output_dir, sample) 

#         if not overwrite:
#             if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
#                 print("File exists "+output_file+", skipping...")
#                 return 0

#         if not os.path.isdir(output_dir):
#             os.makedirs(output_dir)

#         mc_gene_level_worker(allc_file, output_file, genebody=genebody, contexts=contexts)

#         # compress and name them .bgz
#         sp.run("bgzip -f {}".format(output_file), shell=True)
#         sp.run("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)

#     else: 
#         raise ValueError('Invalid convention! choose from ["CEMBA"]!')
    
#     return 0

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
    # parser = create_parser()
    # args = parser.parse_args()

    ti = time.time()

    log = create_logger()

    allc_file = '/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3C_171206/allc/allc_171213_CEMBA_mm_P56_P63_3C_MOp_CEMBA171206_3C_1_CEMBA171206_3C_3_H4_AD002_indexed.tsv.bgz'
    output_file = 'test.tsv'
    bed_file = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/CEMBA_3C_171206/CEMBA_3C_171206_merged.bed'

    mc_region_level_worker(allc_file, output_file, bed_file, contexts=CONTEXTS)
    # mc_region_level_worker(args.input_allc, args.output_file, args.bed_file, contexts=['CH', 'CG'])
    tf = time.time()
    print("time: %s sec" % (tf-ti))


