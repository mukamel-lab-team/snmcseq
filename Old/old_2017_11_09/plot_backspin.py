#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

input_table = './preprocessed/human_hv1_hv2_CGCH_genebody_nmcc_backspin_all.tsv'
tsne_fname = './preprocessed/human_hv1_hv2_CGCH_genebody_nmcc_tsne_all.tsv'
output_plot = './preprocessed/human_hv1_hv2_CGCH_genebody_nmcc_backspin_all.pdf' 

title = 'BackSPIN clustering on genebody mCH-mCG of glial cells'

# make tSNE plot

def plot_tsne_label(df, tx, ty, tc, 
                    output=None, show=True, close=False, 
                    b_ylim=None, title=None, figsize=(8,6)):
    """

    """
    fig, ax = plt.subplots(figsize=figsize)

    for label, df_sub in df.groupby(tc):
        ax.scatter(df_sub[tx], df_sub[ty], s=2, label=label)

    if title:
        ax.set_title(title)

    ax.set_xlabel('tsne_x')
    ax.set_ylabel('tsne_y')
    ax.set_aspect('equal')
    # ax.legend()

    fig.tight_layout()
    if output:
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)

df_cluster = pd.read_table(input_table, header=0)
df_tsne = pd.read_table(tsne_fname, header=0)
df_plot = pd.merge(df_cluster, df_tsne, left_on='Sample', right_on='cells')

plot_tsne_label(df_plot, 'tsne_x', 'tsne_y', 'cluster_ID', 
                title=title,
                output=output_plot, 
                figsize=(6,8))