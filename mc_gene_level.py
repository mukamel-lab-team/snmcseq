#!/usr/bin/python3
#
# Read in allc tables and generate summarized mc levels within gene bodies
#

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
    outdir="./genebody"):
    """
    samples = list of file names of allc files
    genebody = BED file with gene body or other annotation features
    """

    if sample.endswith('_bismark'):
        sample = sample[:-8]

    chromosomes = snmcseq_utils.get_human_chromosomes()

    df_gtf = pd.read_table(genebody)
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[3:])
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    print("mc_gene_level processing: "+sample)

    sample_basename = os.path.basename(sample)
    outfilename = os.path.join(outdir, sample_basename+"_mch_genebody.txt")

    if os.path.isfile(outfilename):
         print("File exists "+outfilename+", skipping...")
         return 0

    outfile_CH = open(outfilename, "w")
    outfile_CH.write("id\tname\tchr\tstart\tend\tstrand\tmc\tc\n")

    for group in df_gtf.groupby(['chr']):
        chrom, df = group
        
        allc = tabix.open(sample+'_bismark/allc_'+sample_basename+'_'+str(chrom)+'.tsv.gz')
        for i, row in df.iterrows():
            records = allc.query(row['chr'], row['start'], row['end'])
            mc, c = snmcseq_utils.tabix_summary(records, context="CH", cap=2)
            outfile_CH.write(row['gene_id'] + "\t" + row['name'] + "\t" + row['chr'] + "\t" + str(row['start']) + "\t" + 
                 str(row['end']) + "\t" + row['strand'] + "\t" + str(mc) + "\t" + str(c) + "\n")

    # for i,row in df_gtf.iterrows():

    #     allc = tabix.open(sample+'_bismark/allc_'+sample_basename+'_'+row['chr']+'.tsv.gz')
    #     # Gene body CH
    #     records = allc.query(row['chr'], row['start'], row['end'])
    #     mc, c = snmcseq_utils.tabix_summary(records, context="CH", cap=2)
    #     outfile_CH.write(row['gene_id'] + "\t" + row['name'] + "\t" + row['chr'] + "\t" + str(row['start']) + "\t" + 
    #        str(row['end']) + "\t" + row['strand'] + "\t" + str(mc) + "\t" + str(c) + "\n")

    return 0

        # procs = min(len(samples), 16)
        # p = mp.Pool(processes=procs)
        # split_samples = np.array_split(samples,procs)
        # pool_results = p.map(process, split_samples)
        # p.close()
        # p.join()


def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", help="sample name")
    parser.add_argument("--genebody", help="file of gene body")
    parser.add_argument("--outdir", help="output directory")
    return parser


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()

    if not args.genebody:
        args.genebody = '/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv'

    if not args.outdir:
        args.outdir = './genebody'

    ti = time.time()
    mc_gene_level(args.sample, genebody=args.genebody, outdir=args.outdir)
    tf = time.time()
    print("time: %s sec" % (tf-ti))


