#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import shutil
import glob
import time
import subprocess as sp
from collections import OrderedDict
import argparse
import logging

import snmcseq_utils
from snmcseq_utils import create_logger
from __init__ import *


# setup ensemble
def initial_setup(ens_path, readme_f, ens_description, 
    meta_fin, meta_fout, allc_paths, genebody_paths, binc_paths, bin_size=BIN_SIZE):
    """
    - Create folder, README.txt
    - Copy over metadata
    - Check if metadata agree with allc files, genebody files, and binc files
    """
    if not os.path.isdir(ens_path):
        # create folder
        os.makedirs(ens_path)
        logging.info("{} created!".format(ens_path))
        
        # create README.txt
        with open(readme_f, 'w') as readme:
            readme.write(ens_description + '\n')
        logging.info("{} created:\n{}".format(readme_f, ens_description[:min(len(ens_description), 50)]))
        
        # copy metadata over
        shutil.copyfile(meta_fin, meta_fout)
        logging.info("{} created from {}".format(meta_fout, meta_fin))
        
        # check if metadata agree with allc/gene_level/binc info
        # for singleton ensemble, metadata samples need to match allc/gene_level/binc
        # for non-singleton ensemble, metadata samples need to be in allc/gene_level/binc
        meta_df = pd.read_table(meta_fin, index_col='Sample')
        cells = sorted(meta_df.index.tolist())
        cells_allc = [os.path.basename(allc_path)[len('allc_'):-len('.tsv.bgz')]
                             for allc_path in allc_paths] 
        cells_genebody = [os.path.basename(genebody_path)[len('genebody_'):-len('.tsv.bgz')]
                             for genebody_path in genebody_paths] 
        cells_binc = [os.path.basename(binc_path)[len('binc_'):-len('_{}.tsv.bgz'.format(bin_size))]
                             for binc_path in binc_paths] 
        
        assert cells == cells_allc
        assert cells == cells_genebody
        assert cells == cells_binc
        logging.info("Mapping summary (cell metadata) matches allc, gene_level, and binc files!")
        
    else:
        # log.info("Error: {} already exists!".format(ens_path))
        raise ValueError("Error: {} already exists!".format(ens_path))

    return cells, cells_allc, cells_genebody, cells_binc 

def pull_genebody_info(ens, ens_genelevel_path, cells_genebody, genebody_paths, 
                contexts=CONTEXTS, to_file=True):
    """
    Pull genebody information from dataset

    each mark (context) per file
    
    return dfs and contexts
    """
    # pull genebody info from dataset
    logging.info("Pulling genebody information ({} cells)...".format(len(cells_genebody)))
    ti = time.time()

    contexts = CONTEXTS
    if not os.path.isdir(ens_genelevel_path):
        # create folder
        os.makedirs(ens_genelevel_path)
        logging.info("{} created!".format(ens_genelevel_path))
        
    output_fnames = [os.path.join(ens_genelevel_path, 'genebody_m{}_{}.tsv'.format(context, ens)) 
                     for context in contexts]

    for i, (cell, genebody_path) in enumerate(zip(cells_genebody, genebody_paths)):
        
        df_gnb = pd.read_table(genebody_path, index_col='gene_id', compression='gzip')
        # df_gnb = pd.read_table(genebody_path, index_col='gene_id', compression='gzip').sort_index()
        
        if i == 0:
            dfs = [pd.DataFrame(index=df_gnb.index)]*len(contexts)
            
        for j, context in enumerate(contexts):
        
            assert dfs[j].index.tolist() == df_gnb.index.tolist()
            
            dfs[j][cell+'_mc'] = df_gnb['m'+context] 
            dfs[j][cell+'_c'] = df_gnb[context] 
            
        logging.info('Loaded cell: {} ({}/{})'.format(cell, i+1, len(cells_genebody)))

    if to_file:
        for df, context, output in zip(dfs, contexts, output_fnames):
            df.to_csv(output, sep='\t', na_rep='NA', index=True, header=True)
            # compress and name them .bgz
            sp.run("bgzip -f {}".format(output), shell=True)
            sp.run("mv {}.gz {}.bgz".format(output, output), shell=True)

            logging.info('Output genebody info to \n{}.bgz'.format(output))

    tf = time.time()
    logging.info("Time spent on pulling genebody information: {} sec".format(tf - ti))

    return dfs, contexts

def pull_binc_info(ens, ens_binc_path, cells_binc, binc_paths, 
                bin_size=BIN_SIZE, contexts=CONTEXTS, to_file=False):
    """
    Pull genebody information from dataset

    each mark (context) per file

    return dfs and contexts

    dfs is indexed by 'chr' and 'bin'
    """
    logging.info("Pulling binc information ({} cells)...".format(len(cells_binc)))
    ti = time.time()

    if not os.path.isdir(ens_binc_path):
        # create folder
        os.makedirs(ens_binc_path)
        logging.info("{} created!".format(ens_binc_path))
        
    output_fnames = [os.path.join(ens_binc_path, 'binc_m{}_{}_{}.tsv'.format(context, ens, bin_size)) 
                     for context in contexts]

    for i, (cell, binc_path) in enumerate(zip(cells_binc, binc_paths)):
        
        df_bin = snmcseq_utils.read_binc(binc_path, compression='gzip')
        
        if i == 0:
            dfs = [pd.DataFrame(index=df_bin.index)]*len(contexts)
            
        for j, context in enumerate(contexts):
        
            assert dfs[j].index.tolist() == df_bin.index.tolist()
            
            dfs[j][cell+'_mc'] = df_bin['m'+context] 
            dfs[j][cell+'_c'] = df_bin[context] 
            
        logging.info('Loaded cell: {} ({}/{})'.format(cell, i+1, len(cells_binc)))

    if to_file:
        for df, context, output in zip(dfs, contexts, output_fnames):
            df.to_csv(output, sep='\t', na_rep='NA', index=True, header=True)
            # compress and name them .bgz
            sp.run("bgzip -f {}".format(output), shell=True)
            sp.run("mv {}.gz {}.bgz".format(output, output), shell=True)

            logging.info('Output binc info to \n{}.bgz'.format(output))

    tf = time.time()
    logging.info("Time spent on pulling binc information: {} sec".format(tf - ti))

    return dfs, contexts


def merge_bins(df, bin_size=10*BIN_SIZE, double_xsize=True, 
               output_file=None):
    """
    Merge bins of BIN_SIZE to n*BIN_SIZE, where n has to be an integer.
    The last incomplete bin for each chromosome is removed.
   
    df has columns: ['chr', 'bin'] (0 based) and [$sample_mc, $sample_c, ....]
    
    return binc file (or choose to save to file)
    """ 
    logging.info("Merging bins to bin_size={}, double_xsize={}".format(bin_size, double_xsize))
    ti = time.time()
    
    chromosomes = snmcseq_utils.get_mouse_chromosomes()
    chrs_all = np.asarray([])
    bins_all = np.asarray([])
    mc_c_all = OrderedDict()
    
    for col in df.columns:
        if col not in ['chr', 'bin']:
            mc_c_all[col] = np.array([])
        
    for chromosome, df_sub in df.groupby('chr'):
        # here -1 is very important!
        bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], bin_size) - 1)

        if double_xsize and chromosome == 'X':
            # here -1 is very important!
            bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], 2*bin_size) - 1)
        
        res = df_sub.groupby(pd.cut(df_sub['bin'], bins)).sum().fillna(0)
        
        chrs = np.asarray([chromosome]*(len(bins)-1))
        bins_all = np.concatenate([bins_all, (bins+1)[:-1]]) # +1 to restore 0-based
        chrs_all = np.concatenate([chrs_all, chrs])
        
        for col in df.columns:
            if col not in ['chr', 'bin']:
                mc_c_all[col] = np.concatenate([mc_c_all[col], res[col]])
        
    # binc
    columns = ['chr', 'bin'] + [key for key in mc_c_all]
    binc = pd.DataFrame(columns=columns)
    binc['chr'] = chrs_all.astype(object)
    binc['bin'] = bins_all.astype(int)
    for key, value in mc_c_all.items():
        binc[key] = value.astype(int) 
    
    if output_file:
        binc.to_csv(output_file, na_rep='NA', sep="\t", header=True, index=False)
        # compress and name them .bgz
        sp.run("bgzip -f {}".format(output_file), shell=True)
        sp.run("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)
        logging.info("Done with binc processing, saving results to: {}.bgz".format(output_file))

    tf = time.time()
    logging.info("Done with binc processing!")
    logging.info("Time spent on pulling binc information: {} sec".format(tf - ti))
    return binc


def main_setup(dataset, ens, ens_description):
    """
    create an singleton ensemble from a dataset:
    - 
    - 
    - 
    """
    # create logger
    log = create_logger()

    log.info("Initiating an singleton ensemble {} from {} \n{}".format(ens, dataset, ens_description))
    tit = time.time()

    ### ---- set up paths ---- ### 
    # path_datasets is the path of "datasets" folder
    bin_size=BIN_SIZE
    path_datasets = PATH_DATASETS
    path_ensembles = PATH_ENSEMBLES

    # dataset_path or dataset_paths is/are the path(s) of specific datasets
    dataset_path = os.path.join(path_datasets, dataset) 
    ens_path = os.path.join(path_ensembles, ens)

    # metadata (mapping summaries)
    meta_fin = os.path.join(dataset_path, 'mapping_summary_{}.tsv'.format(dataset))
    meta_fout = os.path.join(ens_path, 'mapping_summary_{}.tsv'.format(ens))

    # file paths should be sorted!
    # allcs
    allc_paths = sorted(glob.glob(os.path.join(dataset_path, 'allc/allc_*.tsv.bgz')))
    # genebodys 
    genebody_paths = sorted(glob.glob(os.path.join(dataset_path, 'gene_level/genebody_*.tsv.bgz')))
    # bincs
    binc_paths = sorted(glob.glob(os.path.join(dataset_path, 'binc/binc_*_{}.tsv.bgz'.format(bin_size))))

    # README.txt
    readme_f = os.path.join(ens_path, 'README_{}.txt'.format(ens))

    # ensemble genelevel path
    ens_genelevel_path = os.path.join(ens_path, 'gene_level')
    # ensemble binc path
    ens_binc_path = os.path.join(ens_path, 'binc')

    ### ---- end of path set up ---- ### 

    # initial setup
    cells, cells_allc, cells_genebody, cells_binc = initial_setup(
        ens_path, readme_f, ens_description, meta_fin, meta_fout, 
        allc_paths, genebody_paths, binc_paths, bin_size=bin_size)

    # # pull genebody info
    pull_genebody_info(ens, ens_genelevel_path, cells_genebody, genebody_paths, contexts=CONTEXTS)

    # pull binc info
    dfs_binc, contexts = pull_binc_info(ens, ens_binc_path, cells_binc, binc_paths, 
        bin_size=bin_size, contexts=CONTEXTS, to_file=False)

    # reduce binc info to 100kb bins (200kb for chrX)
    for df_binc, context in zip(dfs_binc, contexts):
        df_binc = df_binc.reset_index()
        output_binc = os.path.join(ens_path, 'binc/binc_m{}_{}_{}.tsv'.format(context, ens, 10*bin_size))
        merge_bins(df_binc, bin_size=10*bin_size, double_xsize=True, output_file=output_binc)

    tft = time.time()
    log.info("Ensemble initiation complete: {}".format(ens)) 
    log.info("Total time spent: {} sec".format(tft - tit))


def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", 
        required=True,
        help="Name of the dataset. Eg: CEMBA_3C_171206")
    parser.add_argument("-e", "--ensemble", 
        required=True,
        help="Name of the ensemble. Eg: Ens1")
    parser.add_argument("-m", "--message", 
        required=True,
        nargs="+",
        help="A brief description of the ensemble. Eg: Singleton ensemble of CEMBA_3C_171206 dataset")

    return parser



if __name__ == '__main__':

    # define names and paths 
    parser = create_parser()
    args = parser.parse_args() 

    dataset = args.dataset 
    ens = args.ensemble 
    ens_description = ' '.join(args.message) 

    main_setup(dataset, ens, ens_description)