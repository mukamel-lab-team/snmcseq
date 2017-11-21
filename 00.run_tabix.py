#!/usr/bin/env python3

import os
import glob
import multiprocessing as mp
import numpy as np

import snmcseq_utils 
from snmcseq_utils import create_logger



def make_allc_tabix(sample_path):
    """
    run tabix on all allc files under the path
    args: sample allc file
    """
    if os.path.isdir(sample_path):
        print('Processing: '+sample_path)
        sample = os.path.basename(sample_path)
        for chromosome in snmcseq_utils.get_human_chromosomes() + ['Y']:
            os.system('tabix -f -s 1 -b 2 -e 2 -S 1 ' + 
                    os.path.join(sample_path, 'allc_'+sample+'_'+chromosome+'.tsv.gz')
                    )

            # os.system('bgzip -c allc_'+sample+'_'+chromosome+'.tsv > allc_'+sample+'_'+chromosome+'.tsv.gz')
            # print('allc_'+sample+'_'+chromosome+'.tsv.gz')
            # os.system('tabix -f -s 1 -b 2 -e 2 -S 1 allc_'+sample+'_'+chromosome+'.tsv.gz')
    else:
        print('WARNING: THE BELOW DIRECTORY DOES NOT EXIST.')
        print(sample_path)

    return True

def make_allc_tabixes(allc_dir, nprocs=16):
    """
    Given an allc directory, make tabix files for all samples under the directory
    """

    local_dirs = glob.glob(os.path.join(allc_dir, '*_indexed')) 

    nprocs = min(nprocs, len(local_dirs))
    pool = mp.Pool(processes=nprocs)
    pool_results = [pool.apply_async(make_allc_tabix, 
                                    args=(local_dir, )) 
                    for local_dir in local_dirs]
                    
    pool.close()
    pool.join()

    return pool_results



if __name__ == '__main__':

    # make_allc_tabix('./data/allc/MB_EA/171030_MB_EB_hs_25yr_BA10_FCU2_NA_H9_AD010_indexed')
    # make_allc_tabix('./data/allc/MB_EA/170831_MB_EA_hs_58yr_BA10_pool_2540_AD002_indexed')
    
    ALLC_DIR = './data/allc/MB_EA'

    logger = create_logger()
    logger.info('Begin ...')

    make_allc_tabixes(ALLC_DIR, nprocs=16) 

    logger.info('Done!')
