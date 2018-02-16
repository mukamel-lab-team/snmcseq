"""Config CEMBA scripts
"""
import time
import logging
import glob
import os
import numpy as np
import pandas as pd


# define constant variables
BIN_SIZE = 10000
BIN_SIZE_FEATURE = 10*BIN_SIZE

CONTEXTS = ['CH', 'CG', 'CA']
COMBINED_CONTEXTS_LIST = [['CH', 'CG'], ['CA', 'CG']]


PATH_CEMBA = '/cndd/Public_Datasets/CEMBA/snmCSeq'
PATH_DATASETS = PATH_CEMBA + '/Datasets'
PATH_ENSEMBLES = PATH_CEMBA + '/Ensembles'
PATH_REFERENCES = PATH_CEMBA + '/References'
GENEBODY = PATH_REFERENCES + '/Annotation/gencode.vM16.annotation_genes.tsv'
GENOME_SIZE_FILE = PATH_REFERENCES + '/Genome/mm10.chrom.sizes'

# tSNE
PERPLEXITIES = [20, 30, 40, 50, 100] 
N_PC = 50 
N_DIM = 2

# louvain
K_NN = [5, 10, 15, 20, 30, 50, 100] 

# mysql
DATABASE = 'CEMBA'
CELLS_TABLE_COLS = ['cell_id', 
                     'cell_name', 
                     'dataset', 
                     'cell_type',
                     'global_mCH', 
                     'global_mCG',
                     'global_mCA',
                     'global_mCCC', 
                     'estimated_mCH', 
                     'estimated_mCG',
                     'percent_genome_covered', 
                     'total_reads',
                     'mapped_reads', 
                     'mapping_rate', 
                     'nonclonal_reads', 
                     'percent_nonclonal_rate',
                     'filtered_reads',
                     'filtered_rate',
                     'lambda_mC']

def rename_ms_cols(column_names):
    """
    """
    dict_rename = {'Sample': 'cell_name', 
                   'Total reads': 'total_reads', 
                   'Mapped reads': 'mapped_reads', 
                   'Mapping rate': 'mapping_rate',
                   'Nonclonal reads': 'nonclonal_reads',
                   '% Nonclonal rate': 'percent_nonclonal_rate', 
                   'Filtered reads': 'filtered_reads', 
                   'Filtered rate': 'filtered_rate', 
                   'Lambda mC/C': 'lambda_mC', 
                   'mCCC/CCC': 'global_mCCC', 
                   'mCG/CG': 'global_mCG', 
                   'mCH/CH': 'global_mCH',
                   'Estimated mCG/CG': 'estimated_mCG', 
                   'Estimated mCH/CH': 'estimated_mCH', 
                   '% Genome covered': 'percent_genome_covered'}
    return [dict_rename[col] for col in column_names] 