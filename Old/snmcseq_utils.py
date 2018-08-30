"""
library from Chris's mypy

Fangming edited
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import seaborn as sns

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

def get_id_from_name(gene_name, 
    f_gene_id_to_names='/cndd/fangming/snmcseq_dev/data/references/gene_id_to_names.tsv'):
    """
    get the first matching gene_id from gene_name
    """
    df_gene_id_to_names = pd.read_table(f_gene_id_to_names, index_col='geneID') 
    gene_name = gene_name.upper()
    try:
        gene_id = df_gene_id_to_names[df_gene_id_to_names.geneName==gene_name].index.values[0]
        return gene_id
    except:
        raise ValueError('No matching gene id is found for gene name: %s' % gene_name)


def set_value_by_percentile(this, lo, hi):
    """set `this` below or above percentiles to given values
    this (float)
    lo(float)
    hi(float)
    """
    if this < lo:
        return lo
    elif this > hi:
        return hi
    else:
        return this

def mcc_percentile_norm(mcc, low_p=5, hi_p=95):
    """
    set values above and below specific percentiles to be at the value of percentiles 

    args: mcc, low_p, hi_p  

    return: normalized mcc levels
    """
#   mcc_norm = [np.isnan(mcc) for mcc_i in list(mcc)]
    mcc_norm = np.copy(mcc)
    mcc_norm = mcc_norm[~np.isnan(mcc_norm)]

    lo = np.percentile(mcc_norm, low_p)
    hi = np.percentile(mcc_norm, hi_p)

    mcc_norm = [set_value_by_percentile(mcc_i, lo, hi) for mcc_i in list(mcc)]
    mcc_norm = np.array(mcc_norm)

    return mcc_norm

def plot_tsne_values(df, tx='tsne_x', ty='tsne_y', tc='mCH',
                    s=2,
                    cbar_label=None,
                    output=None, show=True, close=False, 
                    t_xlim=None, t_ylim=None, title=None, figsize=(8,6)):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    """
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.scatter(df[tx], df[ty], s=s, 
        c=mcc_percentile_norm(df[tc].values))
    if title:
        ax.set_title(title)
    ax.set_xlabel(tx)
    ax.set_ylabel(ty)
    ax.set_aspect('auto')
    clb = plt.colorbar(im, ax=ax)
    if cbar_label:
        clb.set_label(cbar_label, rotation=270, labelpad=10)
    if t_xlim:
        ax.set_xlim(t_xlim)
    if t_ylim:
        ax.set_ylim(t_ylim)

    fig.tight_layout()
    if output:
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)



def tsne_and_boxplot(df, tx='tsne_x', ty='tsne_y', tc='mCH', bx='cluster_ID', by='mCH',
                    output=None, show=True, close=False, title=None, figsize=(6,8),
                    t_xlim=None, t_ylim=None, b_xlim=None, b_ylim=None, 
                    low_p=5, hi_p=95):
    """
    boxplot and tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    """
    fig, axs = plt.subplots(2,1,figsize=figsize)

    ax = axs[0]
    im = ax.scatter(df[tx], df[ty], s=2, 
        c=mcc_percentile_norm(df[tc].values, low_p=low_p, hi_p=hi_p))
    # ax.set_xlim([-40, 40])
    # ax.set_ylim([40, 100])
    if title:
        ax.set_title(title)
    ax.set_xlabel('tsne_x')
    ax.set_ylabel('tsne_y')
    ax.set_aspect('auto')
    clb = plt.colorbar(im, ax=ax)
    clb.set_label(tc, rotation=270, labelpad=10)
    if t_ylim:
        ax.set_ylim(t_ylim)
    if t_xlim:
        ax.set_xlim(t_xlim)

    ax = axs[1]
    sns.boxplot(x=bx, y=by, data=df, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if b_ylim:
        ax.set_ylim(b_ylim)
    if b_xlim:
        ax.set_xlim(b_xlim)

    fig.tight_layout()
    if output:
        # output_dir = './preprocessed/marker_genes_%s' % method
        # if not os.path.exists(output_dir):
        #   os.makedirs(output_dir)

        # output_fname = os.path.join(output_dir, '%s_%s.pdf' % (cluster_ID, gene_name))
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)

def myScatter(ax, df, x, y, l, 
              s=20,
              random_state=None,
              legend_mode=0,
              colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'], **kwargs):
    """
    take an axis object and make a scatter plot
    """
    # shuffle (and copy) data
    df = df.sample(frac=1, random_state=random_state)
    # add a color column
    inds, catgs = pd.factorize(df[l])
    df['c'] = [colors[i%len(colors)] if i!=-1 else 'grey' for i in inds]
    # modify label column
    df[l] = df[l].fillna('unlabelled') 
    
    # take care of legend
    for ind, row in df.groupby(l).first().iterrows():
        ax.scatter(row.tsne_x, row.tsne_y, c=row.c, label=ind, s=s)
        
    if legend_mode == 0:
        ax.legend()
    elif legend_mode == -1:
        pass
    elif legend_mode == 1:
        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
              ncol=6, fancybox=False, shadow=False) 
    
    # actual plot
    ax.scatter(df.tsne_x, df.tsne_y, c=df['c'], s=s)
    
    return


def plot_tsne_labels(df, tx='tsne_x', ty='tsne_y', tc='cluster_ID', 
                    legend_mode=0,
                    s=1,
                    random_state=None,
                    output=None, show=True, close=False, 
                    t_xlim=None, t_ylim=None, title=None, figsize=(8,6),
                    legend_loc='lower right',
                    colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9'], **kwargs):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    # avoid gray-like 'C7' in colors
    # color orders are arranged for exci-inhi-glia plot 11/1/2017
    """
    fig, ax = plt.subplots(figsize=figsize)

    myScatter(ax, df, tx, ty, tc,
             s=s,
             random_state=random_state, 
             legend_mode=legend_mode, 
             colors=colors, **kwargs)

    if title:
        ax.set_title(title)
    ax.set_xlabel(tx)
    ax.set_ylabel(ty)
    ax.set_aspect('auto')
    if t_xlim:
        ax.set_xlim(t_xlim)
    if t_ylim:
        ax.set_ylim(t_ylim)

    if output:
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)


# def plot_tsne_labels(df, tx='tsne_x', ty='tsne_y', tc='cluster_ID', 
#                     legend_mode=0,
#                     s=1,
#                     output=None, show=True, close=False, 
#                     t_xlim=None, t_ylim=None, title=None, figsize=(8,6),
#                     legend_loc='lower right',
#                     colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9']):
#     """
#     tSNE plot

#     xlim, ylim is set to facilitate displaying glial clusters only

#     # avoid gray-like 'C7' in colors
#     # color orders are arranged for exci-inhi-glia plot 11/1/2017
#     """
#     fig, ax = plt.subplots(figsize=figsize)

#     for i, (label, df_sub) in enumerate(df.fillna('unlabelled').groupby(tc)):
#         if label == 'unlabelled':
#             ax.scatter(df_sub[tx], df_sub[ty], c='gray',
#                     label=label, s=s) 
#         else:
#             ax.scatter(df_sub[tx], df_sub[ty], c=colors[i%len(colors)], 
#                     label=label, s=s)


#     if title:
#         ax.set_title(title)
#     ax.set_xlabel(tx)
#     ax.set_ylabel(ty)
#     ax.set_aspect('auto')
#     if t_xlim:
#         ax.set_xlim(t_xlim)
#     if t_ylim:
#         ax.set_ylim(t_ylim)

#     if legend_mode == 0:
#         ax.legend()
#         fig.tight_layout()

#     elif legend_mode == -1:
#         fig.tight_layout()
        
#     elif legend_mode == 1:
#         # Shrink current axis's height by 10% on the bottom
#         box = ax.get_position()
#         ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                          box.width, box.height * 0.9])
#         ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
#               ncol=6, fancybox=False, shadow=False) 
#         # add_legend(loc=)

#     if output:
#         fig.savefig(output)
#         print('Saved to ' + output) 
#     if show:
#         plt.show()
#     if close:
#         plt.close(fig)
