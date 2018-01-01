#!/usr/bin/env python3
import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import time
import argparse

import snmcseq_utils


################
#  CREATE THE BINNED DATA FILES FROM ALLC FILES
################


######### INPUT PARAMETERS #################
# species = 'human' # or human
# samples = ['Pool_2256_AD010_indexed_R1'] # LIST OF SAMPLE NAMES
# allc_dir = '/cndd/Public_Datasets/single_cell_methylome/allc_singlecells/hs_MB_EB/' 
# # DIRECTORY WITH ALLC FILES IN IT

# FILE NAME TO ALLC FILES WILL BE CONSTRUCTED AS: allc_dir + "/allc_" + sample[n] + "_" + chromosome + ".tsv"
# OUTPUT FILE NAME: "binc_" + sample[n] + "_" + str(bin_size) + "_" + chromosome + ".tsv",

# species = 'human'
# sample = '/cndd/Public_Datasets/single_cell_methylome/allc_singlecells/hs_MB_EB/Pool_2256_AD006_indexed_R1_bismark'
### ---


### FUNCTION FOR BINNING THE ALLC FILES
def bin_allc(sample_path, 
    bin_size=10000, 
    chromosomes=None, 
    outpath = '.', 
    species='mouse', 
    compressed=False):

    sample = os.path.basename(sample_path)
    if chromosomes == None:
        if species == 'human':
            chromosomes = snmcseq_utils.get_human_chromosomes()
        elif species == 'mouse':
            chromosomes = snmcseq_utils.get_mouse_chromosomes()

    for chromosome in chromosomes:

        fname = sample_path + "/allc_" + sample + "_" + chromosome + ".tsv"

        if compressed:
            os.system("bgzip -cd " + fname + ".gz > " + fname)

        if not os.path.isfile(fname):
            print("bin_allc: " + fname + " does not exist.")
            return

        if not os.path.exists(outpath):
            os.makedirs(outpath)
            
        output_filename = outpath + "/binc_" + sample + "_" + str(bin_size) + "_" + chromosome + ".tsv"
        if os.path.isfile(output_filename):
            print("File exists "+output_filename+", skipping...")
            return 0
        else:
            print("Processing: " + output_filename)

        df = snmcseq_utils.read_allc(fname)
        
        if compressed:
            os.remove(fname)

        if species == 'human':
            bins = np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], bin_size)
        else:
            bins = np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], bin_size)

        # mCG
        df_CG = df.loc[df.context.isin(snmcseq_utils.get_mCG_contexts())]
        groups = df_CG.groupby(pd.cut(df_CG.index, bins))
        mCG = groups.sum().mc.fillna(0)
        CG = groups.sum().c.fillna(0)

        # mCH
        df_CH = df.loc[df.context.isin(snmcseq_utils.get_mCH_contexts())]
        groups = df_CH.groupby(pd.cut(df_CH.index, bins))
        mCH = groups.sum().mc.fillna(0)
        CH = groups.sum().c.fillna(0)

        data = np.array([bins[:len(bins)-1], mCG.values, CG.values, mCH.values, CH.values]).astype(int)
        binned_allc = pd.DataFrame(data.transpose(), columns=['bin','mCG','CG','mCH','CH'])
        binned_allc['chr'] = chromosome
        binned_allc = binned_allc[['chr','bin','mCG','CG','mCH','CH']]

        binned_allc.to_csv(output_filename,
                           na_rep='NA', sep="\t", header=True, index=False)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input directory containing allc files', required=True)
    parser.add_argument('-o', '--output', help='directory storing output files', required=True)
    parser.add_argument('-s', '--species', help='mouse or human', default='mouse')
    parser.add_argument('-bz', '--bin_size', help='bin size', default=10000, type=int)
    parser.add_argument('-c', '--compressed', help='compressed or not', action='store_true') 

    return parser


if __name__ == '__main__':
    
    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()
    bin_allc(args.input, 
            outpath=args.output, 
            species=args.species, 
            bin_size=args.bin_size, 
            compressed=args.compressed)
    tf = time.time()
    print("time: %s sec" % (tf-ti))



