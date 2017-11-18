#!/usr/bin/python3
"""
# Read in allc tables and generate summarized mc levels within gene bodies
"""

import pandas as pd
import tabix
import multiprocessing as mp
import numpy as np
import os
import argparse
import time

import snmcseq_utils

def mc_gene_level(sample,
    genebody='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv',
    outdir="./genebody",
    context='CH'):
    """
    sample = path of allc files of one sample (following certain naming convention)
    genebody = BED file with gene body or other annotation features
    """

    # if sample.endswith('_bismark'):
    #     sample = sample[:-8]

    chromosomes = snmcseq_utils.get_human_chromosomes()

    df_gtf = pd.read_table(genebody)
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[3:])
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    print("mc_gene_level processing: %s %s" % (sample, context))

    sample_basename = os.path.basename(sample)
    outfilename = os.path.join(outdir, sample_basename+("_m%s_genebody.txt" % context))

    if os.path.isfile(outfilename):
         print("File exists "+outfilename+", skipping...")
         return 0

    outfile = open(outfilename, "w")
    outfile.write("id\tname\tchr\tstart\tend\tstrand\tmc\tc\n")

    for group in df_gtf.groupby(['chr']):
        chrom, df = group
        df = df.sort_values('start') 

        allc = tabix.open(sample+'/allc_'+sample_basename+'_'+str(chrom)+'.tsv.gz')
        for i, row in df.iterrows():
            records = allc.query(row['chr'], row['start'], row['end'])
            mc, c = snmcseq_utils.tabix_summary(records, context=context, cap=2)
            outfile.write(row['gene_id'] + "\t" + row['name'] + "\t" + row['chr'] + "\t" + str(row['start']) + "\t" + 
                 str(row['end']) + "\t" + row['strand'] + "\t" + str(mc) + "\t" + str(c) + "\n")

    outfile.close()
    return 0

def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_sample", help="full path of the sample directory of allc files", required=True)
    parser.add_argument("--genebody", help="file of gene body", default='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv')
    parser.add_argument("-o", "--outdir", help="output directory", default='./genebody')
    parser.add_argument("-c", "--context", help="context: CH/CG/...", required=True)
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    ti = time.time()
    mc_gene_level(args.input_sample, genebody=args.genebody, outdir=args.outdir, context=args.context)
    tf = time.time()
    print("time: %s sec" % (tf-ti))


