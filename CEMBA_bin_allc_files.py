#!/usr/bin/env python3
import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import time
import argparse
from collections import OrderedDict
import subprocess as sp

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger


################
#  CREATE THE BINNED DATA FILES FROM ALLC FILES
################



### FUNCTION FOR BINNING THE ALLC FILES
def bin_allc_worker(allc_file, output_file,
    bin_size=BIN_SIZE, 
    contexts=CONTEXTS,
    chromosomes=None, 
    species='mouse',
    compression='gzip'
    ):

    logger = create_logger()
    logger.info("binc processing: {} {}".format(allc_file, contexts))

    if chromosomes == None:
        if species == 'human':
            chromosomes = snmcseq_utils.get_human_chromosomes()
        elif species == 'mouse':
            chromosomes = snmcseq_utils.get_mouse_chromosomes()

    df_allc = snmcseq_utils.read_allc_CEMBA(allc_file, pindex=False, compression=compression)
    df_allc = df_allc.loc[df_allc.chr.isin(chromosomes)]

    chr_allc = np.array([])
    bins_allc = np.array([])
    mc_c_allc = OrderedDict()
    for context in contexts:
        mc_c_allc['m'+context] = np.array([])
        mc_c_allc[context] = np.array([])

    for chromosome, df in df_allc.groupby('chr'):
        # if chromosome == 'X':
        #     if species == 'human':
        #         bins = np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], 2*bin_size)
        #     elif species == 'mouse':
        #         bins = np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], 2*bin_size)
        #     else:
        #         raise ValueError("No such species available: {}".format(species))
        # else:

        # last bin (incomplete) is discarded
        if species == 'human':
            bins = np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], bin_size)
        elif species == 'mouse':
            bins = np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], bin_size)
        else:
            raise ValueError("No such species available: {}".format(species))

        # number of intervals is number of bin points -1 
        chrs = np.asarray([chromosome]*(len(bins)-1))
        bins_allc = np.concatenate([bins_allc, bins[:-1]])
        chr_allc = np.concatenate([chr_allc, chrs])

        # mCG, mCH, mCA
        for context in contexts:
            df_context = df.loc[df.context.isin(snmcseq_utils.get_expanded_context(context))]
            df_mc_c = df_context.groupby(pd.cut(df_context.pos, bins)).sum().fillna(0)[['mc', 'c']]
            mc_c_allc['m'+context] = np.concatenate([mc_c_allc['m'+context], df_mc_c.mc])
            mc_c_allc[context] = np.concatenate([mc_c_allc[context], df_mc_c.c])

    columns = ['chr', 'bin'] + [key for key in mc_c_allc]
    binc = pd.DataFrame(columns=columns)
    binc['chr'] = chr_allc.astype(object)
    binc['bin'] = bins_allc.astype(int)
    for key, value in mc_c_allc.items():
        binc[key] = value.astype(int) 

    binc.to_csv(output_file, na_rep='NA', sep="\t", header=True, index=False)
    logger.info("Done with binc processing: {} {}\nSaving results to: {}".format(allc_file, contexts, output_file))

def bin_allc(allc_file, 
    convention='CEMBA',
    bin_size=BIN_SIZE, 
    contexts=CONTEXTS,
    chromosomes=None, 
    species='mouse',
    compression='gzip',
    overwrite=False
    ):
    """
    set up conventions for output_file
    """

    if convention=='CEMBA':
        CEMBA_DATASETS = PATH_DATASETS 
        allc_file = os.path.abspath(allc_file)
        assert allc_file[:len(CEMBA_DATASETS)] == CEMBA_DATASETS

        dataset, *dis, allc_basename = allc_file[len(CEMBA_DATASETS)+1:].split('/')

        sample = allc_basename[len('allc_'):-len('.tsv.bgz')] 

        output_dir = "{}/{}/binc".format(CEMBA_DATASETS, dataset)
        output_file = "{}/binc_{}_{}.tsv".format(output_dir, sample, bin_size) 

        if not overwrite:
            if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
                logging.info("File exists "+output_file+", skipping...")
                return 0

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        bin_allc_worker(allc_file, output_file,
                bin_size=bin_size, 
                contexts=CONTEXTS,
                chromosomes=chromosomes, 
                species=species,
                compression=compression)

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
    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc", help="allc file path", required=True)
    parser.add_argument('-s', '--species', help='mouse or human', default='mouse')
    parser.add_argument('-bz', '--bin_size', help='bin size', default=BIN_SIZE, type=int)
    parser.add_argument("-c", "--contexts", help="list of contexts: CH/CG/...", nargs='+', default=CONTEXTS)
    parser.add_argument('-cp', '--compression', help='compression type of allc file (bgz is gzip)', default='gzip') 
    parser.add_argument("-chr", "--chromosomes", help="list of chromosomes", nargs='+', default=None)
    parser.add_argument("-f", "--overwrite", 
        action='store_true',
        help="overwrite a file if it exists")
    return parser


if __name__ == '__main__':
    
    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()
    bin_allc(args.input_allc, 
        bin_size=args.bin_size, 
        contexts=args.contexts,
        chromosomes=args.chromosomes, 
        species=args.species,
        compression=args.compression,
        overwrite=args.overwrite
        )
    tf = time.time()
    print("time: %s sec" % (tf-ti))



