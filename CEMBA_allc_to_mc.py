#!/usr/bin/env python3
"""Given allc tables and corresponding DMR bed files, 
prepare mc_single file format (each chromosome) to be ready to be uploaded to mySQL as AnnoJ browser tracks.

Allc -> split allc -> allc2mc -> mc+dmr_to_mc_single
Fangming Xie
"""

import sys
sys.path.insert(0, '/cndd/fangming/CEMBA/snmcseq_dev')

from __init__ import *
from natsort import natsorted

import subprocess as sp
import re
import snmcseq_utils

class_dict = {
    'CGA': 'CG',
    'CGC': 'CG',
    'CGG': 'CG',
    'CGT': 'CG',
    'CGN': 'CG',
    'CAG': 'CHG',
    'CCG': 'CHG',
    'CTG': 'CHG',
    'CAA': 'CHH',
    'CAC': 'CHH',
    'CAT': 'CHH',
    'CCA': 'CHH',
    'CCC': 'CHH',
    'CCT': 'CHH',
    'CTA': 'CHH',
    'CTC': 'CHH',
    'CTT': 'CHH',
} 

def allc_to_mc(allc_table, output_fname=None):
    """Input can be a file path 
    - or an pandas dataframe with columns: ['chr', 'pos', 'strand', 'class', 'mc', 'c']
    
    Return an pandas dataframe with columns: ['assembly', 'position', 'strand', 'class', 'mc', 'c']
    """
    if isinstance(allc_table, str):
        table = snmcseq_utils.read_allc_CEMBA(allc_table, pindex=False, compression='infer')
    else:
        table = allc_table.copy()
        
    table['class'] = table['context'].apply(lambda x: class_dict[x])
    table = table[['chr', 'pos', 'strand', 'class', 'mc', 'c']].rename(columns={
        'chr': 'assembly',
        'pos': 'position',
        'strand': 'strand',
        'class': 'class',
        'mc': 'mc',
        'c': 'h',
    })
    
    if output_fname:
        table.to_csv(output_fname, sep='\t', na_rep='NA', header=True, index=False)
    
    return table

def split_allc(allc_table):
    """Input can be a file path 
    - or an pandas dataframe with columns: ['chr', 'pos', 'strand', 'class', 'mc', 'c']
    """
    if isinstance(allc_table, str):
        table = snmcseq_utils.read_allc_CEMBA(allc_table, pindex=False, compression='infer')
    else:
        table = allc_table.copy()
        
    dict_allc = {}
    for chrom, df_chrom in table.groupby('chr'):
        dict_allc[chrom] = df_chrom
        
    return dict_allc
        
def append_dmr_mc(mc_file, dmr_file, 
                  tmp_file = '/scratch/tmp_dmr_mc'):
    """Find mc in mc_file that are in dmr regions and append them into the mc_files
    Use bash
    """
    if not os.path.isfile(mc_file):
        return
    if not os.path.isfile(dmr_file):
        return
    
    cmd = "cp {} {}".format(mc_file, tmp_file)
    sp.run(cmd, shell=True)
    
    cmd = ("awk \'BEGIN{{FS=\"\\t\"; OFS=\"\\t\"}} NR>1 {{$4=$4\"dmr\"; print $1,$2-1,$2,$0}}\' {2} "
           "| bedtools intersect -wa -u -a - -b {1} "
           "| cut --complement -f 1-3 "
           ">> {0}"
          ).format(mc_file, dmr_file, tmp_file)
    sp.run(cmd, shell=True)
    
    cmd = "rm {}".format(tmp_file)
    sp.run(cmd, shell=True)
    
    return
 
def allc_to_mc_worker(allc_files, output_prefixes):
    """
    Args:
        - allc_files (allc_table 'allc_.tsv' or 'allc_.tsv.gz' or 'allc_.tsv.bgz')
        - output_prefixes (output will be '{}_{}.tsv'.format(output_prefixes, chrom))
    """
    # allc to mc
    for i, (allc_file, output_prefix) in enumerate(zip(allc_files, output_prefixes)):
        logging.info("processing: {}".format(allc_file))
        # split allc
        allc_bychr = split_allc(allc_file)
        for i, (chrom, allc) in enumerate(allc_bychr.items()):
            # allc to mc
            output_file = "{}_{}.tsv".format(output_prefix, chrom)
            allc_to_mc(allc, output_fname=output_file)
    return

def allc_to_mc(allc_files, output_dir):
    """
    Args:
        - allc_files (allc_table 'allc_.tsv' or 'allc_.tsv.gz' or 'allc_.tsv.bgz')
        - output_dir 
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    ### internal
    output_prefixes = [os.path.join(
        output_dir, 
        os.path.basename(allc_file)
        .replace('allc', 'mc')
        .split('.')[0]
        )
        for allc_file in allc_files]
    
    logging.info(np.array(list(zip(allc_files, output_prefixes))))
    
    # work
    allc_to_mc_worker(allc_files, output_prefixes)

    return


if __name__ == '__main__':

    log = snmcseq_utils.create_logger()

    allc_files = natsorted(glob.glob(os.path.join(
        './allc_multimodal_v2_*.tsv'))) 

    output_dir = '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/mc_tracks'

    allc_dmr_to_mc_single(allc_files, dmr_files, output_dir, use_dmrs)



    
