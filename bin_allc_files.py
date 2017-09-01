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


def get_mCH_contexts():
    contexts = []
    for base1 in ['A','C','T']:
        for base2 in ['A','C','G','T']:
            contexts.append('C' + base1 + base2)
    return contexts

def get_mCG_contexts():
    return ['CGA','CGC','CGG','CGT']

def get_mouse_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,20)]
    if include_x:
        chromosomes.append('X')
    return chromosomes

def get_human_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,23)]
    if include_x:
        chromosomes.append('X')
    return chromosomes

def get_sample_from_path(path):
    sample = os.path.basename(path)
    if sample.endswith('_bismark'):
        sample = sample[:-8]
    return sample

def get_chrom_lengths_mouse():
    return {'1': 195471971, '2': 182113224, '3': 160039680, '4': 156508116, '5': 151834684, 
            '6': 149736546, '7': 145441459, '8': 129401213, '9': 124595110, '10': 130694993, 
            '11': 122082543, '12': 120129022, '13': 120421639, '14': 124902244, '15': 104043685, 
            '16': 98207768, '17': 94987271, '18': 90702639, '19': 61431566, 'X': 171031299, 'Y': 91744698}

def get_chrom_lengths_human():
    return {'1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555, '5': 181538259, '6': 170805979, 
            '7': 159345973, '8': 145138636, '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309, 
            '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345, '17': 83257441, '18': 80373285, 
            '19': 58617616, '20': 64444167, '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415}

### FUNCTION FOR BINNING THE ALLC FILES
def bin_allc(sample, path='.', bin_size=10000, chromosomes=None, outpath = '.', species='mouse', compressed=False):

    if chromosomes == None:
        if species == 'human':
            chromosomes = get_human_chromosomes()
        else:
            chromosomes = get_mouse_chromosomes()

    for chromosome in chromosomes:

        fname = path + "/allc_" + sample + "_" + chromosome + ".tsv"

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

        df = snmcseq_utils.read_allc(fname)
        
        if compressed:
            os.remove(fname)

        # for bin_size in bin_sizes:

        if species == 'human':
            bins = np.arange(0, get_chrom_lengths_human()[chromosome], bin_size)
        else:
            bins = np.arange(0, get_chrom_lengths_mouse()[chromosome], bin_size)

        # mCG
        df_CG = df.loc[df.context.isin(get_mCG_contexts())]
        groups = df_CG.groupby(pd.cut(df_CG.index, bins))
        mCG = groups.sum().mc.fillna(0)
        CG = groups.sum().c.fillna(0)

        # mCH
        df_CH = df.loc[df.context.isin(get_mCH_contexts())]
        groups = df_CH.groupby(pd.cut(df_CH.index, bins))
        mCH = groups.sum().mc.fillna(0)
        CH = groups.sum().c.fillna(0)

        data = np.array([bins[:len(bins)-1], mCG.values, CG.values, mCH.values, CH.values]).astype(int)
        binned_allc = pd.DataFrame(data.transpose(), columns=['bin','mCG','CG','mCH','CH'])
        binned_allc['chr'] = chromosome
        binned_allc = binned_allc[['chr','bin','mCG','CG','mCH','CH']]

        binned_allc.to_csv(output_filename,
                           sep="\t", header=False, index=False)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input directory containing allc files', required=True)
    parser.add_argument('-o', '--output', help='directory storing output files', required=True)
    parser.add_argument('-s', '--species', help='mouse or human', default='mouse')
    parser.add_argument('-bz', '--bin_size', help='bin size', default=10000, type=int)
    parser.add_argument('-c', '--compressed', help='compressed or not', action='store_true') 

    return parser


if __name__ == '__main__':
    
    # if not os.path.exists(sample):
    #     os.makedirs(sample)


    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()
    bin_allc(get_sample_from_path(args.input), path=args.input, 
            outpath=args.output, species=args.species, bin_size=args.bin_size, compressed=args.compressed)
    tf = time.time()
    print("time: %s sec" % (tf-ti))


### PROCESS THE WILL CALL THE BINNING FUNCTION IN A PARALLEL FASHION
# def process(samples):
#     for sample in samples:
#         print(sample)
#         if not os.path.exists(sample):
#             os.makedirs(sample)
#         bin_allc(sample, path=allc_dir+sample, 
#                 outpath=sample, species=species, bin_size=100000, compressed=False)
#     return True


# procs = min(16, len(samples))

# p = mp.Pool(processes=procs)
# split_samples = np.array_split(samples,procs)
# pool_results = p.map(process, split_samples)
# p.close()
# p.join()

