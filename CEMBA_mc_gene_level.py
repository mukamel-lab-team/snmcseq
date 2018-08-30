#!/usr/bin/python3
"""
# Read in allc tables (with tabix open) and generate summarized mc levels within gene bodies
"""

import pandas as pd
import tabix
import numpy as np
import os
import argparse
import time
import subprocess as sp

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger

def mc_gene_level_worker(allc_file, output_file,
    genebody=GENEBODY,
    contexts=CONTEXTS):
    """
    allc_file 
    genebody anntotation file:
    """
    logger = create_logger()

    chromosomes = snmcseq_utils.get_chromosomes(SPECIES)

    df_gtf = pd.read_table(genebody)
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[len('chr'):])
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    logger.info("mc_gene_level processing: {} {}".format(allc_file, contexts))


    outfile = open(output_file, "w")
    columns = ['gene_id']
    for context in contexts:
        columns += ['m{}'.format(context), context]

    outfile.write('\t'.join(columns)+'\n')

    allc = tabix.open(allc_file)
    for i, row in df_gtf.iterrows():
        row_out = [row.gene_id]
        records = list(allc.query(row['chr'], row['start'], row['end']))
        for context in contexts:
            mc, c = snmcseq_utils.tabix_summary(records, context=context, cap=2) # remove sites with total_c > 2
            row_out += [str(mc), str(c)] 

        outfile.write('\t'.join(row_out)+'\n')

    outfile.close()
    logger.info("Done with mc_gene_level processing: {} {}\n Saving results to {}".format(allc_file, contexts, output_file))

    return 0


def mc_gene_level(allc_file,
    contexts=CONTEXTS, 
    genebody=GENEBODY,
    convention='CEMBA', 
    overwrite=False):

    """
    set up conventions for output_file
    """

    if convention=='CEMBA':
        CEMBA_DATASETS = PATH_DATASETS 
        allc_file = os.path.abspath(allc_file)
        assert allc_file[:len(CEMBA_DATASETS)] == CEMBA_DATASETS

        dataset, *dis, allc_basename = allc_file[len(CEMBA_DATASETS)+1:].split('/')


        sample = allc_basename[len('allc_'):-len('.tsv.bgz')] 

        output_dir = "{}/{}/gene_level".format(CEMBA_DATASETS, dataset)
        output_file = "{}/genebody_{}.tsv".format(output_dir, sample) 

        if not overwrite:
            if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
                print("File exists "+output_file+", skipping...")
                return 0

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        mc_gene_level_worker(allc_file, output_file, genebody=genebody, contexts=contexts)

        # compress and name them .bgz
        try:
            sp.run("bgzip -f {}".format(output_file), shell=True)
            sp.run("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)
        except:
            sp.call("bgzip -f {}".format(output_file), shell=True)
            sp.call("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)

    else: 
        raise ValueError('Invalid convention! choose from ["CEMBA"]!')
    
    return 0

def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc", help="allc file path", required=True)
    parser.add_argument("-c", "--contexts", help="list of contexts: CH/CG/...", nargs='+', default=CONTEXTS)
    parser.add_argument("-g", "--genebody", help="file of gene body", default=GENEBODY)
    parser.add_argument("-f", "--overwrite", 
        action='store_true',
        help="overwrite a file if it exists")
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()
    mc_gene_level(args.input_allc, contexts=args.contexts, genebody=args.genebody, overwrite=args.overwrite)
    tf = time.time()
    print("time: %s sec" % (tf-ti))


