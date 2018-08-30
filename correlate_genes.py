"""
Calculates the top 50 most highly correlated genes for each gene using spearman correlation.

Input: Gene x Cells table containing mc and c columns for each cell.
Output: N*50 x 3 .tsv file containing gene1, gene2, and correlation value in each row. N = number of genes. 

"""

# data i/o
import time
import csv
import argparse

# the usual
import numpy as np
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')


from __init__ import *
import snmcseq_utils


def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-ei", "--ensemble_id",
                        required=True,
                        type=int,
                        help="Ensemble ID, Eg: 1, 2 ,3, 4, etc...")
    parser.add_argument("-p", "--pearson",
                        action="store_true",
                        help="Toggle ranking of rows for pearson correlation. Else, spearman correlation is calculated.")
    parser.add_argument("-c", "--context",
                        default='CH',
                        help="Methylation context (CG, CH, CA). Default = CH.")
    parser.add_argument("-n", "--ngenes",
                        default="50",
                        type=int,
                        help="Number of top correlating genes to keep. Default = 50.")
    return parser


def correlation(ENSEMBLE_ID, SPEARMAN, CONTEXT, NUM_GENES):
    # ## PREPROCESS DATA

    start_time = time.time()
    print("Reading Data", flush=True)
    ensemble = "Ens" + ENSEMBLE_ID
    in_file = "/cndd/tishihar/single_cell_methylome/correlations/gene_level/genebody_m{0}_{1}.tsv".format(CONTEXT, ensemble)
    gene_x_cells = pd.read_table(in_file)

    print("Creating mcc dataframe..", flush=True)
    df_mc = gene_x_cells.filter(regex="_mc$", axis=1)
    df_c = gene_x_cells.filter(regex="_c$", axis=1)
    df_mc.columns = df_mc.columns + 'c'
    df_c.columns = df_mc.columns
    df_mcc = df_mc/df_c
    df_mcc.insert(0, 'gene_id', gene_x_cells['gene_id'])


    #### Impute NaN values with mean for each gene ####

    print("Imputing data..", flush=True)
    df_mcc = df_mcc.loc[df_mcc.loc[:, df_mcc.columns!='gene_id'].count(axis=1) > 0]
    df_mcc.reset_index(inplace=True, drop=True)
    means = df_mcc.loc[:, df_mcc.columns!='gene_id'].mean(axis=1)
    fill_value = pd.DataFrame({col: means for col in df_mcc.columns if col != 'gene_id'})
    df_mcc.fillna(fill_value, inplace=True)


    # Rank each row for spearman correlation and save binary file
    correlation_type = "pearson"
    if !PEARSON:
        print("Ranking each row for Spearman correlation..", flush=True)
        df_mcc = df_mcc.iloc[:, df_mcc.columns != 'id'].rank(axis=1)
        correlation_type = "spearman"


    # ## CORRELATION MATRIX

    print("Building correlation matrix..", flush=True)
    corr_matrix = np.corrcoef(df_mcc.as_matrix())
    corr_matrix[np.isnan(corr_matrix)] = 0 # Get rid of NaN values

    print("Sorting correlation matrix..", flush=True)

    # Argsort and get top 50 for each row (descending order)
    corr_matrix_argsorted = (-corr_matrix).argsort(axis=1)
    corr_matrix_argsorted_top_genes = corr_matrix_argsorted[:, :NUM_GENES]

    # Sort and get top 50 for each row (descending order)
    corr_matrix_sorted = np.sort(corr_matrix, axis=1)[:, ::-1]
    corr_matrix_sorted_top_genes = corr_matrix_sorted[:, :NUM_GENES]


    # Match gene_ids to argsorted correlation matrix.
    print("Building final table", flush=True)
    corr_matrix_sorted_gene_ids_flat = np.array([gene_x_cells['gene_id'].to_dict().get(i,-1) for i in range(corr_matrix_argsorted_top_genes.min(), corr_matrix_argsorted_top_genes.max()+1)])
    corr_matrix_sorted_gene_ids = corr_matrix_sorted_gene_ids_flat[(corr_matrix_argsorted_top_genes - corr_matrix_argsorted_top_genes.min())]

    # Zip the gene_id and correlation value for each gene. e.g. gene1: [(gene2_0, corr_0), (gene2_1, corr_1) ... (gene2_51149, corr_51149)]
    genes_and_correlation_dict = {gene_id: list(zip(corr_matrix_sorted_gene_ids[i], corr_matrix_sorted_top_genes[i])) for i, gene_id in enumerate(corr_matrix_sorted_gene_ids[:,0])}

    # Create the intended format.
    final_list = []
    for gene1, correlations in genes_and_correlation_dict.items():
        for correlation in correlations:
            final_list.append((gene1, correlation[0], correlation[1]))

    output_file = "/cndd/tishihar/single_cell_methylome/correlations/{0}_{1}_correlated_{2}_genes_m{3}.tsv".format(ensemble, correlation_type, NUM_GENES, CONTEXT)

    with open(output_file, 'w') as f:
        csv_out = csv.writer(f, delimiter='\t')
        csv_out.writerow(['gene1', 'gene2', 'correlation'])
        for row in final_list:
            csv_out.writerow([row[0], row[1], round(row[2], 6)])

    print("DONE! File saved at {}".format(output_file), flush=True)
    print("Total time (minutes): " + str(round((time.time()-start_time)/60, 2)), flush=True)

    return


if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()

    ENSEMBLE_ID = str(args.ensemble_id)
    PEARSON = args.pearson
    CONTEXT = args.context
    NUM_GENES = args.ngenes

    correlation(ENSEMBLE_ID, SPEARMAN, CONTEXT, NUM_GENES)
