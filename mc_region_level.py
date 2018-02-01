#!/usr/bin/python3
"""
# Read in allc tables and generate summarized mc levels within regions specified by bed file 
"""

import pandas as pd
import tabix
import multiprocessing as mp
import numpy as np
import os
import argparse
import glob
import time
from collections import OrderedDict

import snmcseq_utils

def mc_region_level_worker(allc_file_dict,
    bed_file,
    output_file,
    contexts=['CH']):

    """
    allc_files is a dictionary of {'chrom': allc_chrom_file_name, ...} for one sample 
    """
    chromosomes = snmcseq_utils.get_human_chromosomes()

    df_gtf = pd.read_table(bed_file, header=None)
    df_gtf.columns = ['chr', 'start', 'end'] + [i+1 for i in range(df_gtf.shape[1]-3)]
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[len('chr'):]) # remove chr
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)] # keep those that are in the chromosomes only

    print("mc_region_level processing: {} {}".format(allc_file_dict, contexts))


    if os.path.isfile(output_file):
         print("File exists "+output_file+", skipping...")
         return 0

    outfile = open(output_file, "w")

    columns = ['chr', 'start', 'end']
    for context in contexts:
        columns += ['m{}'.format(context), context]

    outfile.write('\t'.join(columns)+'\n')

    for chrom, df in df_gtf.groupby(['chr']):
        df = df.sort_values('start') 
        allc_file = allc_file_dict[chrom]
        print(chrom)
        allc = tabix.open(allc_file)
        for i, row in df.iterrows():
            row_out = [str(row.chr), str(row.start), str(row.end)]
            records = list(allc.query(row['chr'], row['start'], row['end']))
            for context in contexts:
                mc, c = snmcseq_utils.tabix_summary(records, context=context, cap=0) # remove sites with total_c > 2
                row_out += [str(mc), str(c)] 

            outfile.write('\t'.join(row_out)+'\n')

    outfile.close()
    return 0


def mc_region_level(allc_dir,
    bed_file,
    output_file,
    contexts=['CH'],
    convention='classical'):
    """
    given allc_dir (containing just one sample), create allc_file_dict

    """
    if convention == 'classical':
        allc_files = glob.glob('{}/*.tsv.gz'.format(allc_dir))
        allc_file_dict = OrderedDict({f.split('.')[0].split('_')[-1]: f for f in allc_files})

        mc_region_level_worker(allc_file_dict, bed_file, output_file, contexts=contexts)

        return

def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc_dir", help="Path of the allc directory that contain allc tables (from different chromosomes) for one sample", required=True)
    parser.add_argument("-b", "--bed_file", help="--bed_file", required=True)
    parser.add_argument("-o", "--output_file", help="output file", required=True)
    parser.add_argument("-c", "--contexts", nargs="+", help="list of contexts: CH/CG/...", default=['CH', 'CG'])
    return parser

if __name__ == '__main__':
    # parser = create_parser()
    # args = parser.parse_args()

    ti = time.time()
    parser = create_parser()
    args = parser.parse_args()

    allc_dir = args.input_allc_dir
    bed_file = args.bed_file 
    output_file = args.output_file 
    contexts = args.contexts

    mc_region_level(allc_dir, bed_file, output_file, contexts=contexts)

    tf = time.time()
    print("time: %s sec" % (tf-ti))


