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

def mc_region_level(allc_file_dict,
    bed_file,
    output_file,
    contexts=['CH']):

    """
    allc_files is a dictionary of {'chorom': allc_chrom_file_name, ...}
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

# def create_parser():
#     """

#     """
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-i", "--input_sample", help="full path of the sample directory of allc files", required=True)
#     parser.add_argument("--genebody", help="file of gene body", default='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv')
#     parser.add_argument("-o", "--outdir", help="output directory", default='./genebody')
#     parser.add_argument("-c", "--context", help="context: CH/CG/...", required=True)
#     return parser


if __name__ == '__main__':
    # parser = create_parser()
    # args = parser.parse_args()

    ti = time.time()

    allc_dir = '/cndd/Public_Datasets/single_cell_methylome/allc_combined/human'
    # hL2/3 - 13, hPv1 - 4

    allc_files = glob.glob(os.path.join(allc_dir, 'allc_human_cluster13_*.tsv.gz'))
    allc_files = [os.path.join(allc_dir, allc_file) for allc_file in allc_files]
    allc_file_dict = OrderedDict({f.split('.')[0].split('_')[-1]: f for f in allc_files})
    # print(allc_file_dict)

    bed_file = '/cndd/junhao/genomes/hg19/genomicRegions/UCSC_repeatMasker_LINE.bed'
    output_file = './hL2-3_LINE.tsv'
    contexts = ['CH', 'CG']

    mc_region_level(allc_file_dict, bed_file, output_file, contexts=contexts)

    tf = time.time()
    print("time: %s sec" % (tf-ti))


