"""
library from Chris's mypy

Fangming edited
"""

import pandas as pd
import logging

def create_logger(name='log'):
    """
    args: logger name

    return: a logger object
    """
    logging.basicConfig(
        format='%(asctime)s %(message)s', 
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    return logging.getLogger(name)



def get_mCH_contexts(wildcard=True):
    contexts = []
    for base1 in ['A','C','T']:
        for base2 in ['A','C','G','T']:
            contexts.append('C' + base1 + base2)
    if wildcard:
        return contexts+['CAN', 'CTN', 'CCN']
    else:
        return contexts 

def get_mCG_contexts(wildcard=True):
    if wildcard:
        return ['CGA','CGC','CGG','CGT'] + ['CGN']
    else:
        return ['CGA','CGC','CGG','CGT']


def get_mouse_chromosomes(include_x=True, include_chr=False):
    chromosomes = [str(x) for x in range(1,20)]
    if include_x:
        chromosomes.append('X')
    if not include_chr:
        return chromosomes
    else:
        return ['chr'+chrom for chrom in chromosomes]

def get_human_chromosomes(include_x=True, include_chr=False):
    chromosomes = [str(x) for x in range(1,23)]
    if include_x:
        chromosomes.append('X')
    if not include_chr:
        return chromosomes
    else:
        return ['chr'+chrom for chrom in chromosomes]

def get_chrom_lengths_mouse():
    return {'1': 195471971, '2': 182113224, '3': 160039680, '4': 156508116, '5': 151834684, 
            '6': 149736546, '7': 145441459, '8': 129401213, '9': 124595110, '10': 130694993, 
            '11': 122082543, '12': 120129022, '13': 120421639, '14': 124902244, '15': 104043685, 
            '16': 98207768, '17': 94987271, '18': 90702639, '19': 61431566, 'X': 171031299, 'Y': 91744698}

# hg38
# def get_chrom_lengths_human():
#     return {'1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555, '5': 181538259, '6': 170805979, 
#             '7': 159345973, '8': 145138636, '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309, 
#             '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345, '17': 83257441, '18': 80373285, 
#             '19': 58617616, '20': 64444167, '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415}

# hg19
def get_chrom_lengths_human(
    genome_size_fname='/cndd/fangming/iGenome/hg19/hg19.chrom.sizes'):  
    """
    """
    srs_gsize = pd.read_table(genome_size_fname, header=None, index_col=0, squeeze=True)
    srs_gsize = srs_gsize.loc[get_human_chromosomes(include_chr=True)]
    # remove leading 'chr'
    srs_gsize.index = [idx[len('chr'):] for idx in srs_gsize.index]
    return srs_gsize

def tabix_summary(records, context="CH", cap=0):

    mc = 0
    c = 0

    if context == "CH":
    	# 15 contexts 
        contexts = get_mCH_contexts()
    elif context == "CG":
    	# 5 contexts 
        contexts = get_mCG_contexts()
    elif context == "CA":
        contexts = ["CAA","CAC","CAG","CAT"]
    elif context == "CT":
        contexts = ["CTA","CTC","CTG","CTT"]
    elif context == "CAG":
        contexts = ["CAG"]
    elif context == "CAC":
        contexts = ["CAC"]
    else:
        raise ValueError('Invalid context.')

    if cap > 0:
        for record in records:
            if record[3] in contexts:
                if int(record[5]) <= cap:
                    mc += int(record[4])
                    c += int(record[5])
    else:

        for record in records:
            if record[3] in contexts:
                mc += int(record[4])
                c += int(record[5])

    return mc, c


def read_allc(fname, position_as_index=True, compressed=False):
    if compressed:
        os.system("bgzip -cd " + fname + ".gz > " + fname)

    if position_as_index == True:
        df = pd.read_csv(fname, sep="\t", index_col=1, skiprows=1,
                         names=['chr','pos','strand','context','mc','c','methylated'])
    else:
        df = pd.read_csv(fname, sep="\t", skiprows=1,
                         names=['chr','pos','strand','context','mc','c','methylated'])

    if compressed:
        os.remove(fname)
    return df

# def read_gencode_human(version='v19', pc=False):
# 	# pc = protein coding
#     prefix = '/cndd/projects/Public_Datasets/references/hg19/transcriptome/'
#     if pc:
#         fname= prefix+'gencode.'+version+'.annotation_genes_pc_mypy.tsv'
#     else:
#         fname= prefix+'gencode.'+version+'.annotation_genes_mypy.tsv'
#     return pd.read_csv(fname, sep="\t")
