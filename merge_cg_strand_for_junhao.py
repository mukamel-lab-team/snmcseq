#!/usr/bin/env python3

from __init__ import *
from natsort import natsorted

import snmcseq_utils


def merge_cg_strand(df):
    """Merge allCG table to allCG (one strand) table
    allCG table: chr, pos, strand, context, mc, c, methylation
    """
    # df = df[df['context'].isin(['CGA', 'CGG', 'CGT', 'CGC', 'CGN'])].copy()

    df['pos'] -= (df['strand'].values == '-').astype(int)

    res = df[['chr', 'pos', 'mc', 'c']].groupby(['chr', 'pos']).sum()
    res['strand'] = '+'
    res['context'] = 'CGN'
    res['methylation'] = 1
    res = res[['strand', 'context', 'mc', 'c', 'methylation']].reset_index()
    return res

if __name__ == '__main__':

    log = snmcseq_utils.create_logger()


    files = [
            '/cndd/junhao/Dnmt3aKO/methylome_FC_P0_P13_P39/DMV/allc_tmp.txt',
            ]
    # files

    for i, file in enumerate(files):
        logging.info("Processing file: {} ({}/{})".format(file, i+1, len(files)))
        df = pd.read_table(file, header=None, 
                        names=['chr', 'pos', 'strand', 'context', 'mc', 'c', 'methylation'], 
                        dtype={'chr': str})

        res = merge_cg_strand(df)
        res.to_csv(os.path.splitext(file)[0]+'_mergestrands.tsv', sep='\t', na_rep='NA', 
                    header=False, index=False)
        # res.to_csv('/cndd/fangming/test_mergestrands.tsv', sep='\t', na_rep='NA', 
        #             header=False, index=False)
