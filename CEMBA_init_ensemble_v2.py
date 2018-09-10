#!/usr/bin/env python3

# import matplotlib
# matplotlib.use('Agg')

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

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger
from snmcseq_utils import chunks 
from snmcseq_utils import compress 
import CEMBA_preprocess_bins 
import CEMBA_run_tsne
import CEMBA_clustering_louvain_jaccard
import CEMBA_autoannotate
import CEMBA_update_mysql
import CEMBA_marker_genes
import CEMBA_marker_genes_mysql


# setup ensemble
def initial_setup(ens_path, readme_f, ens_name, ens_datasets, ens_description, 
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
            readme.write("Ensemble name: " + ens_name +'\n')
            readme.write("Ensemble datasets: {}".format(ens_datasets)+'\n')
            readme.write("Total number of cells: {}".format(len(allc_paths))+'\n')

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

# setup ensemble
def initial_setup_nonsingleton(ens_path, readme_f, ens_name, ens_datasets, ens_description, ens_cells, 
    meta_fins, meta_fout, allc_paths, genebody_paths, binc_paths, bin_size=BIN_SIZE):
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
            readme.write("Ensemble name: " + ens_name +'\n')
            readme.write("Ensemble (distinct) datasets: {}".format(ens_datasets)+'\n')
            readme.write("Total number of cells: {}".format(len(ens_cells))+'\n')


        logging.info("{} created:\n{}".format(readme_f, ens_description[:min(len(ens_description), 50)]))
        
        # copy metadata over
        meta_dfs = [pd.read_table(meta_fin, index_col='Sample') for meta_fin in meta_fins]
        meta_df = pd.concat(meta_dfs)
        meta_df = meta_df[meta_df.index.isin(ens_cells)]
        meta_df.to_csv(meta_fout, sep="\t", na_rep="NA", index=True, header=True)
        logging.info("{} created from {}".format(meta_fout, meta_fins))
        
        # check if metadata agree with allc/gene_level/binc info
        # metadata samples need to match allc/gene_level/binc
        cells = meta_df.index.tolist()
        cells_allc = [os.path.basename(allc_path)[len('allc_'):-len('.tsv.bgz')]
                             for allc_path in allc_paths]
        cells_genebody = [os.path.basename(genebody_path)[len('genebody_'):-len('.tsv.bgz')]
                             for genebody_path in genebody_paths] 
        cells_binc = [os.path.basename(binc_path)[len('binc_'):-len('_{}.tsv.bgz'.format(bin_size))]
                             for binc_path in binc_paths] 
        
        assert set(ens_cells) == set(cells)
        assert set(ens_cells) == set(cells_allc)
        assert set(ens_cells) == set(cells_genebody)
        assert set(ens_cells) == set(cells_binc)
        logging.info("Mapping summary (cell metadata) matches allc, gene_level, and binc files!")
        
    else:
        # log.info("Error: {} already exists!".format(ens_path))
        raise ValueError("Error: {} already exists!".format(ens_path))

    return cells, cells_allc, cells_genebody, cells_binc 


def pull_genebody_info(ens, ens_genelevel_path, cells_genebody, genebody_paths, 
                contexts=CONTEXTS, to_file=False):
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
            dfs = []
            for context in contexts:
                dfs.append(pd.DataFrame(index=df_gnb.index))
            
        for j, context in enumerate(contexts):
        
            assert dfs[j].index.tolist() == df_gnb.index.tolist()
            
            dfs[j][cell+'_mc'] = df_gnb['m'+context] 
            dfs[j][cell+'_c'] = df_gnb[context] 
            
        logging.info('Loaded cell: {} ({}/{})'.format(cell, i+1, len(cells_genebody)))

    if to_file:
        for df, context, output in zip(dfs, contexts, output_fnames):
            df.to_csv(output, sep='\t', na_rep='NA', index=True, header=True)
            # compress and name them .bgz
            try:
                sp.run("bgzip -f {}".format(output), shell=True)
                sp.run("mv {}.gz {}.bgz".format(output, output), shell=True)
            except:
                sp.call("bgzip -f {}".format(output), shell=True)
                sp.call("mv {}.gz {}.bgz".format(output, output), shell=True)
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
        
    output_fnames = [os.path.join(ens_binc_path, 'binc_m{}_{}_{}.tsv'.format(context, bin_size, ens)) 
                     for context in contexts]

    for i, (cell, binc_path) in enumerate(zip(cells_binc, binc_paths)):
        
        df_bin = snmcseq_utils.read_binc(binc_path, compression='gzip')
        
        if i == 0:
            dfs = []
            for context in contexts:
                dfs.append(pd.DataFrame(index=df_bin.index))
            # dfs = [pd.DataFrame(index=df_bin.index)]*len(contexts) # never use this 
        # if i < 10: # useful for debug 
        for j, context in enumerate(contexts):
        
            assert dfs[j].index.tolist() == df_bin.index.tolist()
            
            dfs[j][cell+'_mc'] = df_bin['m'+context] 
            dfs[j][cell+'_c'] = df_bin[context] 
            
        logging.info('Loaded cell: {} ({}/{})'.format(cell, i+1, len(cells_binc)))

    if to_file:
        for df, context, output in zip(dfs, contexts, output_fnames):
            df.to_csv(output, sep='\t', na_rep='NA', index=True, header=True)
            # compress and name them .bgz
            try:
                sp.run("bgzip -f {}".format(output), shell=True)
                sp.run("mv {}.gz {}.bgz".format(output, output), shell=True)
            except:
                sp.call("bgzip -f {}".format(output), shell=True)
                sp.call("mv {}.gz {}.bgz".format(output, output), shell=True)

            logging.info('Output binc info to \n{}.bgz'.format(output))

    tf = time.time()
    logging.info("Time spent on pulling binc information: {} sec".format(tf - ti))

    return dfs, contexts


def merge_bins(df, bin_size=10*BIN_SIZE, double_xsize=True, species=SPECIES, 
               output_file=None):
    """
    Merge bins of BIN_SIZE to n*BIN_SIZE, where n has to be an integer.
    The last incomplete bin for each chromosome is removed.
   
    df has columns: ['chr', 'bin'] (0 based) and [$sample_mc, $sample_c, ....]
    
    return binc file (or choose to save to file)
    """ 
    logging.info("Merging bins to bin_size={}, double_xsize={}".format(bin_size, double_xsize))
    ti = time.time()
    
    chromosomes = snmcseq_utils.get_chromosomes(species)
    chrs_all = np.asarray([])
    bins_all = np.asarray([])
    mc_c_all = OrderedDict()
    
    for col in df.columns:
        if col not in ['chr', 'bin']:
            mc_c_all[col] = np.array([])
        
    for chromosome, df_sub in df.groupby('chr'):
        # here -1 is very important!
        if species == 'mouse':
            bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], bin_size) - 1)
        elif species == 'human': 
            bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], bin_size) - 1)
        else:
            raise ValueError("species has to be mouse or human")

        if double_xsize and chromosome == 'X':
            # here -1 is very important!
            if species == 'mouse':
                bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], 2*bin_size) - 1)
            elif species == 'human': 
                bins = (np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], 2*bin_size) - 1)
            else:
                raise ValueError("species has to be mouse or human")

        
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

    # set chr bin as index
    binc = binc.set_index(['chr', 'bin'])
    
    if output_file:
        binc.to_csv(output_file, na_rep='NA', sep="\t", header=True, index=True)
        # compress and name them .bgz
        try:
            sp.run("bgzip -f {}".format(output_file), shell=True)
            sp.run("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)
        except:
            sp.call("bgzip -f {}".format(output_file), shell=True)
            sp.call("mv {}.gz {}.bgz".format(output_file, output_file), shell=True)

        logging.info("Done with binc processing, saving results to: {}.bgz".format(output_file))

    tf = time.time()
    logging.info("Done with binc processing!")
    logging.info("Time spent on pulling binc information: {} sec".format(tf - ti))
    return binc


def main_setup(dataset, ens, ens_description, ens_name=None, ens_datasets=None, species=SPECIES):
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
    logging.info("Setting up paths and names...")
    if not ens_datasets:
        ens_datasets = [dataset] # for singleton ensemble
    if not ens_name:
        ens_name = dataset # for singleton ensemble

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
        ens_path, readme_f, ens_name, ens_datasets, ens_description, meta_fin, meta_fout, 
        allc_paths, genebody_paths, binc_paths, bin_size=bin_size)

    # # # pull genebody info
    # pull_genebody_info(ens, ens_genelevel_path, cells_genebody, genebody_paths, contexts=CONTEXTS)

    # pull binc info
    logging.info("Pulling binc info...")
    dfs_binc, contexts = pull_binc_info(ens, ens_binc_path, cells_binc, binc_paths, 
        bin_size=bin_size, contexts=CONTEXTS, to_file=False)

    # reduce binc info to 100kb bins (200kb for chrX)
    logging.info("Merging binc info into {} bins...".format(BIN_SIZE_FEATURE))
    dfs_merged_binc = []
    for df_binc, context in zip(dfs_binc, contexts):
        df_binc = df_binc.reset_index()
        # added 18-04-04
        if snmcseq_utils.isrs2(dataset):
            df_binc = df_binc[df_binc['chr'].isin(snmcseq_utils.get_chromosomes(species, include_x=False))]
        # added 18-04-04
        output_binc = os.path.join(ens_path, 'binc/binc_m{}_{}_{}.tsv'.format(context, BIN_SIZE_FEATURE, ens))
        merged_binc = merge_bins(df_binc, bin_size=10*bin_size, double_xsize=True, output_file=output_binc)
        dfs_merged_binc.append(merged_binc)
    del dfs_binc

    # preprocess context
    logging.info("Preprocessing binc info for tSNE and clustering...")
    dfs_nmcc = []
    for df_merged_binc, context in zip(dfs_merged_binc, contexts):
        df_nmcc = CEMBA_preprocess_bins.preproc_bins(ens, bin_size=BIN_SIZE_FEATURE, context=context, df=df_merged_binc, to_file=True)
        dfs_nmcc.append(df_nmcc)
    del dfs_merged_binc

    # preprocess combined contexts 
    combined_contexts_list = COMBINED_CONTEXTS_LIST
    dfs_nmcc_combined = []
    for comb_contexts in combined_contexts_list:
        dfs_temp = []
        for context in comb_contexts:
            dfs_temp.append(dfs_nmcc[contexts.index(context)]) 
        df_nmcc_combined = CEMBA_preprocess_bins.preproc_bins_combine_contexts(ens, dfs=dfs_temp, bin_size=BIN_SIZE_FEATURE, contexts=comb_contexts, to_file=True)
        dfs_nmcc_combined.append(df_nmcc_combined)
    del dfs_nmcc
    del dfs_nmcc_combined

    '''In the future, if memory space is not enough, we could choose not giving results to (dfs_binc, dfs_merged_binc, dfs_nmcc, dfs_nmcc_combined),
        and read them from files instead to save memory space. 
    '''
    
    # use (dfs_nmcc, contexts) and (dfs_nmcc_combined, combined_contexts_list)
    # tSNE use from nmcc files
    logging.info("Running tSNE...")
    CEMBA_run_tsne.run_tsne_CEMBA(ens, perps=PERPLEXITIES, n_pc=N_PC, n_dim=N_DIM)
    # louvain clustering
    logging.info("Running clustering...")
    # added 18-04-04
    ks = []
    for k in K_NN:
        if k < (len(cells)/2):
            ks.append(k)
    # added 18-04-04
    CEMBA_clustering_louvain_jaccard.run_louvain_CEMBA(ens, ks=ks, n_pc=N_PC)
    # annotation
    logging.info("Running cluster annotation...")
    CEMBA_autoannotate.run_autoannotate_CEMBA(ens)

    # update to mysql
    logging.info("Upload info to mySQL database...")
    CEMBA_update_mysql.upload_ensemble_level(ens, ens_name, ens_datasets, database=DATABASE)

    # added 18-07-10
    # need ensemble level mysql set up
    logging.info("Running cluster marker genes...")
    CEMBA_marker_genes.find_marker_genes_CEMBA(ens, 
        context='CH', clsts='auto', p_putative=0.10)
    # upload mysql 
    CEMBA_marker_genes_mysql.upload_marker_genes(ens, context='CH', database=DATABASE)

    tft = time.time()
    log.info("Ensemble initiation complete: {}".format(ens)) 
    log.info("Total time spent: {} sec".format(tft - tit))


def main_setup_nonsingleton(ens, ens_name, ens_description, ens_datasets=None, ens_sql=None, ens_cells=None, species=SPECIES):
    """
    There should be 1 not None value for (ens_datasets, ens_sql, ens_cells) 
    """

    # create logger
    log = create_logger()

    # use the option parsing in to generate ensemble. (list of cell names and distinct datasets)
    # get ens_datasets and ens_cells
    database = DATABASE
    engine = CEMBA_update_mysql.connect_sql(database)
    if ens_datasets:
        ens_sql = '''SELECT cell_name FROM cells WHERE dataset IN {}'''.format(tuple(ens_datasets))
        ens_cells = pd.read_sql(ens_sql, engine)['cell_name'].tolist()
    elif ens_sql:
        ens_cells = pd.read_sql(ens_sql, engine)['cell_name'].tolist()
        sql = '''SELECT DISTINCT(dataset) FROM cells WHERE cell_name IN {}'''.format(tuple(ens_cells))
        ens_datasets = pd.read_sql(sql, engine)['dataset'].tolist()
    elif ens_cells:
        sql = '''SELECT DISTINCT(dataset) FROM cells WHERE cell_name IN {}'''.format(tuple(ens_cells))
        ens_datasets = pd.read_sql(sql, engine)['dataset'].tolist() 
    else:
        raise ValueError("No valid option chosen to create an ensemble")
    ens_datasets = sorted(ens_datasets)
    ens_cells = sorted(ens_cells)
    sql = '''SELECT cell_name, dataset FROM cells WHERE cell_name IN {}'''.format(tuple(ens_cells))
    df_cells = pd.read_sql(sql, engine)
    # create ensemble with ens, ens_name, ens_description, ens_datasets, and ens_cells
    log.info("Initiating an non-singleton ensemble {} ({})\nDistinctive datasets: {} \n{}".format(ens, ens_name, ens_datasets, ens_description))
    tit = time.time()

    ### ---- set up paths ---- ### 
    logging.info("Setting up paths and names...")
    # path_datasets is the path of "datasets" folder
    bin_size=BIN_SIZE
    path_datasets = PATH_DATASETS
    path_ensembles = PATH_ENSEMBLES

    # dataset_path or dataset_paths is/are the path(s) of specific datasets
    dataset_paths = [os.path.join(path_datasets, dataset) for dataset in ens_datasets] 
    ens_path = os.path.join(path_ensembles, ens)

    # metadata (mapping summaries)
    meta_fins = [os.path.join(dataset_path, 'mapping_summary_{}.tsv'.format(dataset)) 
                for (dataset, dataset_path) in zip(ens_datasets, dataset_paths)]
    meta_fout = os.path.join(ens_path, 'mapping_summary_{}.tsv'.format(ens))

    # file paths should be sorted!
    allc_paths = []
    genebody_paths = []
    binc_paths = []
    for dataset, df_cells_sub in df_cells.groupby('dataset'):
        dataset_path = os.path.join(path_datasets, dataset)
        for cell_name in df_cells_sub['cell_name']:
            # allcs
            allc_paths.append(os.path.join(dataset_path, 'allc/allc_{}.tsv.bgz'.format(cell_name)))
            # genebodys 
            genebody_paths.append(os.path.join(dataset_path, 'gene_level/genebody_{}.tsv.bgz'.format(cell_name)))
            # bincs
            binc_paths.append(os.path.join(dataset_path, 'binc/binc_{}_{}.tsv.bgz'.format(cell_name, bin_size)))

    # README.txt
    readme_f = os.path.join(ens_path, 'README_{}.txt'.format(ens))

    # ensemble genelevel path
    ens_genelevel_path = os.path.join(ens_path, 'gene_level')
    # ensemble binc path
    ens_binc_path = os.path.join(ens_path, 'binc')

    ### ---- end of path set up ---- ### 

    # initial setup
    cells, cells_allc, cells_genebody, cells_binc = initial_setup_nonsingleton(
        ens_path, readme_f, ens_name, ens_datasets, ens_description, ens_cells, 
        meta_fins, meta_fout, 
        allc_paths, genebody_paths, binc_paths, bin_size=bin_size)

    # # # pull genebody info
    # pull_genebody_info(ens, ens_genelevel_path, cells_genebody, genebody_paths, contexts=CONTEXTS)

    # pull binc info
    n_chunk = 1000 
    logging.info("BEGIN PULLING BINC INFO... ({} at a time, {} in total)".format(n_chunk, len(binc_paths)))
    dfs_binc = []
    for i, (cells_binc_chunk, binc_paths_chunk) in enumerate(zip(chunks(cells_binc, n_chunk), chunks(binc_paths, n_chunk))):
        logging.info("PULLING PROGRESS: {}-{}/{}".format(i*n_chunk+1, (i+1)*n_chunk, len(binc_paths)))
        # pull one chunk
        dfs_binc_chunk, contexts = pull_binc_info(ens, ens_binc_path, 
            cells_binc_chunk, binc_paths_chunk, 
            bin_size=bin_size, contexts=CONTEXTS, to_file=False)

        # reduce binc info to 100kb bins (200kb for chrX)
        for i, (df_binc, context) in enumerate(zip(dfs_binc_chunk, contexts)):
            dfs_binc_chunk[i] = merge_bins(dfs_binc_chunk[i].reset_index(), bin_size=BIN_SIZE_FEATURE, double_xsize=True, output_file=None)

        # store (reduced chunk) in dfs_binc 
        dfs_binc.append(dfs_binc_chunk)

    # combine dfs_binc and save to files (1 context per file) 
    logging.info("BEGIN COMBINING BINC CHUNKS...")
    dfs_merged_binc = []
    for i, context in enumerate(contexts):
        dfs_tmp = [dfs_binc_chunk[i] for dfs_binc_chunk in dfs_binc]
        df_binc = pd.concat(dfs_tmp, axis=1)
        # added 18-04-04
        include_RS2 = False
        for dataset in ens_datasets:
            if snmcseq_utils.isrs2(dataset):
                include_RS2 = True
        if include_RS2:
            df_binc = df_binc.reset_index()
            df_binc = df_binc[df_binc['chr'].isin(snmcseq_utils.get_chromosomes(species, include_x=False))]
            df_binc = df_binc.set_index(['chr', 'bin'])
        # added 18-04-04
        # use for further analysis
        dfs_merged_binc.append(df_binc)
        # save output
        output_binc = os.path.join(ens_path, 'binc/binc_m{}_{}_{}.tsv'.format(context, BIN_SIZE_FEATURE, ens))
        df_binc.to_csv(output_binc, sep="\t", na_rep="NA", header=True, index=True)
        # compress and name them .bgz
        compress(output_binc)
        logging.info("Done with binc processing, saving results to: {}.bgz".format(output_binc))
    del dfs_binc  
    del dfs_tmp
    del df_binc

    # preprocess bins 
    for i, (df_merged_binc, context) in enumerate(zip(dfs_merged_binc, contexts)):
        logging.info("Preprocessing binc info for tSNE and clustering...")
        CEMBA_preprocess_bins.preproc_bins(ens, bin_size=BIN_SIZE_FEATURE, context=context, df=df_merged_binc, to_file=True)
    del dfs_merged_binc 

    # preprocess combined contexts 
    combined_contexts_list = COMBINED_CONTEXTS_LIST
    for comb_contexts in combined_contexts_list:
        dfs_tmp = []
        for context in comb_contexts:
            nmcc_file = os.path.join(ens_path, 'binc/binc_m{}_{}_nmcc_{}.tsv'.format(context, BIN_SIZE_FEATURE, ens))
            df_tmp = pd.read_table(nmcc_file, index_col=['chr', 'bin'], dtype={'chr': object, 'bin': int})
            dfs_tmp.append(df_tmp) 
        CEMBA_preprocess_bins.preproc_bins_combine_contexts(ens, dfs=dfs_tmp, bin_size=BIN_SIZE_FEATURE, contexts=comb_contexts, to_file=True)

    # use (dfs_nmcc, contexts) and (dfs_nmcc_combined, combined_contexts_list)

    # tSNE use from nmcc files
    logging.info("Running tSNE...")
    CEMBA_run_tsne.run_tsne_CEMBA(ens, perps=PERPLEXITIES, n_pc=N_PC, n_dim=N_DIM)
    # louvain clustering
    logging.info("Running clustering...")
    # added 18-04-04
    ks = []
    for k in K_NN:
        if k < (len(cells)/2):
            ks.append(k)
    # added 18-04-04
    CEMBA_clustering_louvain_jaccard.run_louvain_CEMBA(ens, ks=ks, n_pc=N_PC)
    # annotation
    logging.info("Running cluster annotation...")
    CEMBA_autoannotate.run_autoannotate_CEMBA(ens)

    # update to mysql
    logging.info("Upload info to mySQL database...")
    CEMBA_update_mysql.upload_ensemble_level(ens, ens_name, ens_datasets, database=DATABASE)

    # added 18-07-10
    # need ensemble level mysql set up
    logging.info("Running cluster marker genes...")
    CEMBA_marker_genes.find_marker_genes_CEMBA(ens, 
        context='CH', clsts='auto', p_putative=0.10)
    # upload mysql 
    CEMBA_marker_genes_mysql.upload_marker_genes(ens, context='CH', database=DATABASE)

    tft = time.time()
    log.info("Ensemble initiation complete: {}".format(ens)) 
    log.info("Total time spent: {} sec".format(tft - tit))
    return


def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--singleton", 
        action='store_true',
        help="whether or not this is a singleton ensemble. If true, need dataset, ensemble, and message arguments")
    # required by both
    parser.add_argument("-ei", "--ensemble_id", 
        required=True,
        type=int,
        help="Ensemble ID. Eg: 1, 2, 3, 4, etc...")
    parser.add_argument("-m", "--message", 
        required=True,
        nargs="+",
        help="A brief description of the ensemble. Eg: Singleton ensemble of CEMBA_3C_171206 dataset")
    # required by singleton ensemble
    parser.add_argument("-d", "--dataset", 
        help="Name of the dataset. Eg: CEMBA_3C_171206")
    # required by non-singleton ensemble
    parser.add_argument("-en", "--ensemble_name", 
        default=None,
        help="Name of the ensemble. Eg: CEMBA_3C")
    # one of the following required by non-singleton ensemble
    parser.add_argument("--ensemble_datasets", 
        nargs='+',
        help="List of datasets included in the ensemble. Eg: CEMBA_3C_171206 CEMBA_3C_171207")
    parser.add_argument("--ensemble_cells", 
        nargs='+',
        help="List of cell names included in the ensemble (space separated). Eg: 171213_CEMBA_mm_P56_P63_3C_MOp_CEMBA171207_3C_3_CEMBA171207_3C_4_E3_AD010_indexed...")
    parser.add_argument("--ensemble_sql",
        nargs='+', 
        help="A sql query string that specify the cells Eg: SELECT cell_name FROM cells WHERE dataset LIKE 'CEMBA_4B_%'")

    return parser

if __name__ == '__main__':

    # define names and paths 
    parser = create_parser()
    args = parser.parse_args() 

    if args.singleton:
        """
        singleton ensemble
        """
        if not args.dataset:
            parser.error("--singleton requires --dataset")

        dataset = args.dataset 
        ens = "Ens{}".format(args.ensemble_id) 
        ens_description = ' '.join(args.message) 
        ens_name = args.ensemble_name

        main_setup(dataset, ens, ens_description, ens_name=ens_name)

    else:
        """
        combined ensemble
        """
        if not args.ensemble_name:
            parser.error("Non-singleton ensemble requires --ensemble_name")
        if args.dataset:
            parser.error("Non-singleton ensemble doesn't need --dataset, (which is for singleton ensemble).")

        ens = "Ens{}".format(args.ensemble_id) 
        ens_description = ' '.join(args.message) 
        ens_name = args.ensemble_name

        ens_datasets = args.ensemble_datasets 
        ens_cells = args.ensemble_cells
        ens_sql = args.ensemble_sql
        if ens_sql:
            ens_sql = ' '.join(ens_sql)
            ens_sql = ens_sql.replace("%", "%%")

        options = (isinstance(ens_datasets, list), 
                    isinstance(ens_sql, str),
                    isinstance(ens_cells, list)) 

        if sum(options) != 1:
            raise ValueError("Please choose 1 and 1 only from (ensemble_datasets, ensemble_sql, ensemble_cells) to define ensemble")

        kwargs = {}
        if ens_datasets:
            kwargs['ens_datasets'] = ens_datasets
        elif ens_sql:
            kwargs['ens_sql'] = ens_sql 
        elif ens_cells:
            kwargs['ens_cells'] = ens_cells
        assert len(kwargs) == 1

        # parsing 1 option into main_setup ...
        main_setup_nonsingleton(ens, ens_name, ens_description, **kwargs)



    
