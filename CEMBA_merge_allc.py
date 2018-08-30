#!/usr/bin/env python3

import queue 
# import logging
import multiprocessing as mp
import argparse

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger
from CEMBA_update_mysql import connect_sql


def encode_allc_chrom(chrom):
    """give every chromosome an integer name to facilitate sorting (as CEMBA order)
    """
    trans_dict={'L': -4,
                'M': -3, 
                'X': -2, 
                'Y': -1, 
                }
    try:
        chrom = int(chrom)
    except:
        chrom = trans_dict[chrom]
    return chrom

def decode_allc_chrom(chrom_code):
    """give every chromosome code back to chrom (as CEMBA order)
    """
    if chrom_code == -4:
        chrom = 'L' 
    elif chrom_code == -3:
        chrom = 'M' 
    elif chrom_code == -2:
        chrom = 'X' 
    elif chrom_code == -1:
        chrom = 'Y' 
    else:
        chrom = int(chrom_code)
    return chrom

def merge_allc_in_chunks(allc_paths, context='CG', chunksize=100000, n_chunk1=50, n_chunk2=20):
    """Merge allc tables given the allc_files
    Allc files are assumed to have CEMBA format (no header, all chromosomes in one file, and bgzipped)
    """
    iter_allcs = [snmcseq_utils.read_allc_CEMBA(allc_path, chunksize=chunksize, pindex=False)
               for allc_path in allc_paths]
    n_allcs = len(allc_paths)
    merged_cks = []
    i = 0
    ti = time.time()
    while True: # iterate over all iterators once 
                # end while loop if empty for all of them 
        i += 1
        if i%5 == 0:
            logging.info("({}) Merging progress (outer): {}*{}, {} seconds".format(n_allcs, i, chunksize, time.time()-ti))
            ti = time.time()

        # establish queue by iterating over all file once
        q = queue.Queue()

        # read phase  (50 at a time)
        j = 0
        tj = time.time()
        for iter_allc_ck in snmcseq_utils.chunks(iter_allcs, n_chunk1):
            j += 1
            if j%5 == 0:
                logging.info("({}) Merging progress (inner): {}*{}, {} seconds".format(n_allcs, j, n_chunk1, time.time()-tj))
                tj = time.time()
            
            dfs_ck = []
            for iter_allc in iter_allc_ck: # iter_allc could be empty, but j still counts up
                try:
                    df_ck = next(iter_allc)
                    df_ck = df_ck.loc[df_ck.context.isin(snmcseq_utils.get_expanded_context(context)), 
                                      ['chr', 'pos', 'mc', 'c']]
                    dfs_ck.append(df_ck)
                    empty_in = False
                except:
                    pass
                
            if not dfs_ck:
                pass
            else: # concat and merge phase
                df_ck = pd.concat(dfs_ck)
                df_ck['chr_code'] = df_ck['chr'].apply(encode_allc_chrom)
                df_ck = df_ck.set_index(['chr_code', 'pos'])[['mc', 'c']]
                # merge 
                merged_ck_tmp = df_ck.groupby(['chr_code', 'pos']).sum()
                # enqueue
                q.put(merged_ck_tmp)
         
        logging.info("({}) Merging progress (inner): {}*{}, {} seconds".format(n_allcs, j, n_chunk1, time.time()-tj))
            
        # dequeue
        if q.empty(): # end the while-loop
            logging.info("({}) Merging progress (outer): {}*{}, {} seconds".format(n_allcs, i, chunksize, time.time()-ti))
            break  
        else: # concat and merge phase (second merge)
            merged_ck = queue_merge(q, n_chunk2)
            # enqueue
            merged_cks.append(merged_ck)
            
    return merged_cks 


def queue_merge(q, n_chunk):
    """Queue merge
    
    Arguments: q (a queue object with dataframes)
    """
    i = 0
    # get n_chunk out if not empty
    dfs = [q.get() for i in range(n_chunk) if not q.empty()]
    ti = time.time()
    while not q.empty():
        i += 1
        if i%10 == 0:
            logging.info("Queue merging progress: {}*{}, {} seconds".format(i, n_chunk, time.time()-ti))
            ti = time.time()

        # merge them and put back in queue
        df = pd.concat(dfs).groupby(['chr_code', 'pos']).sum()
        q.put(df)

        # get n_chunk out if not empty
        dfs = [q.get() for i in range(n_chunk) if not q.empty()]

    # merge them 
    df_final = pd.concat(dfs).groupby(['chr_code', 'pos']).sum()

    logging.info("Queue merging progress: {}*{}, {} seconds".format(i, n_chunk, time.time()-ti))
    
    return df_final

def merge_allc(allc_paths, output_fname, 
	context='CG', chunksize=100000, n_chunk1=50, n_chunk2=20, compression=True):
	"""
	merge allc tables with CEMBA convention (no header)

	phase1: merge in chunks 
		- merge a chunksize from each allc_files into 1 file in 2 steps:
			- read and merge 50 files * chunksize -> one
			- queue merge the resulting chunk -> one
		- get merged_cks as results 

	phase2: queque merge merged chunks
	"""
	ti = time.time()
	logging.info("Merging allc files (# of files: {}, context: {}, output: {})".format(len(allc_paths), context, output_fname))

	# phase 1
	logging.info("({}) Merging phase 1...".format(len(allc_paths)))
	merged_cks = merge_allc_in_chunks(allc_paths, 
		context=context, chunksize=chunksize, n_chunk1=n_chunk1, n_chunk2=n_chunk2)

	# establish queue
	q = queue.Queue()
	for mck in merged_cks:
	    q.put(mck)
	del merged_cks
	    
	logging.info("({}) Number of dataframes after Phase 1: {}".format(len(allc_paths), q.qsize()))

	# phase 2
	logging.info("({}) Merging phase 2...".format(len(allc_paths)))
	df_final = queue_merge(q, n_chunk2)

	# organize results after phase 2
	logging.info("({}) Organizing results...".format(len(allc_paths)))
	# read allcg
	path_allcg = os.path.join(PATH_REFERENCES, 'Genome/mm10_all_cg.tsv')
	allcg = pd.read_table(path_allcg, dtype={'chr': object})
	allcg['chr_code'] = allcg['chr'].apply(encode_allc_chrom)
	# merge on chr and pos
	df_final = df_final.reset_index()
	df_final = pd.merge(df_final, allcg, on=['chr_code', 'pos'], how='left')
	# select columns
	df_final = df_final.loc[df_final['chr'].isin(snmcseq_utils.get_mouse_chromosomes()+['Y', 'M']), 
	             ['chr', 'pos', 'strand', 'context', 'mc', 'c']]
	df_final['methylated'] = 1

	# saving results
	logging.info("({}) Output shape: {}".format(len(allc_paths), df_final.shape))
	# save to file
	logging.info("({}) Saving results to {}".format(len(allc_paths), output_fname))
	df_final.to_csv(output_fname, sep='\t', header=False, index=False, na_rep='NA')

	# added 07/20/2018
	if compression:
		snmcseq_utils.compress(output_fname, suffix='gz')

	logging.info("({}) Total time spent: {} seconds".format(len(allc_paths), time.time()-ti))
	return


def merge_allc_CEMBA(ens, context='CG', cluster_type='cluster_mCHmCG_lv_npc50_k30', database=DATABASE, nprocs=2,
	chunksize=1000000, n_chunk1=50, n_chunk2=20):
	"""Merge allc tables from the same cluster for a specific ensemble
	Arguments: 
		- ens: ensemble name, eg: Ens7
		- cluster_type: the column name in mySQL ensemble table	
		- database: CEMBA 
	returns: 
		- save merged allc files to PATH_ENSEMBLES/$ens/allc_merged
	"""

	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	engine = connect_sql(database)
	# create folder 
	directory = os.path.join(ens_path, 'allc_merged')
	if not os.path.isdir(directory):
		os.makedirs(directory)
		logging.info("Created dir: {}".format(directory))

	# get a clustering result
	annot_type = 'annotation_' + cluster_type[len('cluster_'):]
	sql = '''SELECT cell_name, dataset, {}, {} 
			FROM cells 
			RIGHT JOIN {} 
			ON cells.cell_id = {}.cell_id'''.format(cluster_type, annot_type, ens, ens)
	df_cluster = pd.read_sql(sql, engine, index_col='cell_name')
	df_cluster.columns = ['dataset', 'cluster', 'annotation'] 
	output_summary = os.path.join(ens_path, 'allc_merged/summary_allc_merged_m{}_{}_{}.tsv'.format(context, cluster_type, ens))
	df_cluster.to_csv(output_summary, sep='\t', na_rep='NA', header=True, index=True)
	logging.info("Saved summary file to {}".format(output_summary))

	n_clusters = len(np.unique(df_cluster['cluster'].values))
	# parallelize 
	# more info
	# group allcs for each cluster_id
	allc_paths_all = []
	cluster_id_all = []
	output_fname_all = []
	# get all information
	for cluster_id, df_sub in df_cluster.groupby('cluster'): 
		allc_paths = [os.path.join(PATH_DATASETS, '{}/allc/allc_{}.tsv.bgz').format(dataset, cell) 
                      for (dataset, cell) in zip(df_sub.dataset, df_sub.index)]
		output_fname = os.path.join(ens_path, 'allc_merged/allc_merged_m{}_{}_{}_{}.tsv'.format(context, cluster_type, cluster_id, ens))

		cluster_id_all.append(cluster_id)
		output_fname_all.append(output_fname)
		allc_paths_all.append(allc_paths)


	# set up parallel merging 
	nprocs = min(nprocs, n_clusters)
	logging.info("""Begin merging allc files:{}, {}\n
				Number of processes:{}\n
				Number of clusters to merge:{}\n
				""".format(ens, cluster_type, nprocs, n_clusters))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(merge_allc, 
									args=(allc_paths, output_fname), 
									kwds={'context': context, 
										'chunksize': chunksize,
										'n_chunk1': n_chunk1,
										'n_chunk2': n_chunk2,
										}) 
					for allc_paths, output_fname in zip(allc_paths_all, output_fname_all)]
					
	pool.close()
	pool.join()

	logging.info('Done!')

	return pool_results

def merge_allc_parallel(allc_paths_all, output_fname_all, context='CG', nprocs=2,
	chunksize=1000000, n_chunk1=50, n_chunk2=20):
	"""Merge allc tables in parallel 
	Arguments: 
		- allc_paths_all: a list of list of allc paths
		- output_fname_all: a list of output files
		- nprocs: number of processes
	returns: 
		- save merged allc files to PATH_ENSEMBLES/$ens/allc_merged
	"""

	n_clusters = len(allc_paths_all)
	assert n_clusters == len(output_fname_all) 

	# set up parallel merging 
	nprocs = min(nprocs, n_clusters)
	logging.info("""Begin merging allc files. 
				Number of processes:{}
				Number of output files to merge:{}
				""".format(nprocs, n_clusters))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(merge_allc, 
									args=(allc_paths, output_fname), 
									kwds={'context': context, 
										'chunksize': chunksize,
										'n_chunk1': n_chunk1,
										'n_chunk2': n_chunk2,
										}) 
					for allc_paths, output_fname in zip(allc_paths_all, output_fname_all)]
					
	pool.close()
	pool.join()

	logging.info('Done!')

	return pool_results



def create_parser():
	"""
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-ei", "--ensemble_id", 
		required=True,
		type=int,
		help="Ensemble ID. Eg: 1, 2, 3, 4, etc...")
	parser.add_argument("-cl", "--cluster_type", 
		default='cluster_mCHmCG_lv_npc50_k30', 
		help="Cluster type shown in the CEMBA database ensemble table. eg: cluster_mCHmCG_lv_npc50_k30")
	parser.add_argument("-c", "--context", 
		default='CG',
		help="context (CG, CH, CA or...)")
	parser.add_argument("-n", "--nprocs", 
		default=2, 
		type=int,
		help="number of proceses to use")
	return parser

if __name__ == '__main__':

	log = create_logger()

	parser = create_parser()
	args = parser.parse_args() 

	ens = 'Ens{}'.format(args.ensemble_id)

	merge_allc_CEMBA(ens, context=args.context, 
		cluster_type=args.cluster_type, database=DATABASE, nprocs=args.nprocs)

