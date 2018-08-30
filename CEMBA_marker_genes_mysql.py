#!/usr/bin/env python3

from __init__ import *
import sqlalchemy as sa
from natsort import natsorted 

import snmcseq_utils
import CEMBA_update_mysql


def define_marker_gene_table(metadata, ens):
    """
    """
    from sqlalchemy.dialects.mysql import TINYINT
    
    # 3 column primary key
    table = sa.Table('{}_cluster_marker_genes'.format(ens), metadata, 
                    sa.Column('clustering', sa.String(255), primary_key=True),
                    sa.Column('cluster', TINYINT, primary_key=True),
                    sa.Column('rank', TINYINT, primary_key=True),
                    sa.Column('gene_id', sa.String(255)),
                    )
    return table


def upload_marker_genes(ens, context, database=DATABASE):
    """
    """
    logging.info("Upload marker genes: {} {}".format(ens, context))

    table_name = '{}_cluster_marker_genes'.format(ens)
    engine = CEMBA_update_mysql.connect_sql(database)
    f_gene = os.path.join(PATH_REFERENCES, 'Annotation', 'gencode.vM16.annotation_genes.tsv')
    
    # create the mysql table
    metadata = sa.MetaData(engine)
    marker_gene_table = define_marker_gene_table(metadata, ens)
    metadata.create_all()

    # get table from files
    files = natsorted(glob.glob(os.path.join(PATH_ENSEMBLES, ens, 'cluster', '*.{}_marker_genes').format(context)))

    # figure out corresponding clustering in mysql
    clst_types = []
    for file in files:
        tmp = os.path.basename(file).split('_')
        clst_type = '_'.join([tmp[5]] + tmp[1:4])
        clst_types.append(clst_type)

    # gene_name to id
    df_genes = pd.read_table(f_gene, index_col='gene_name')
    df_genes = df_genes[~df_genes.index.duplicated(keep='first')]

    # clustering, cluster, rank, gene_id
    dfs = []
    for i, (file, clst_type) in enumerate(zip(files, clst_types)):
        df = pd.read_table(file)
        df['rank'] = df.index.values + 1
        df = df.melt(id_vars=['rank'], var_name='cluster', value_name='gene_name')
        df['clustering'] = clst_type
        df['gene_id'] = df['gene_name'].apply(lambda x: df_genes.loc[x, 'gene_id'])
        df = df[['clustering', 'cluster', 'rank', 'gene_id']]
        dfs.append(df)

    dfs = pd.concat(dfs, ignore_index=True) # ignore index!!!

    # insert into mysql table
    CEMBA_update_mysql.insert_into(engine, table_name, dfs, verbose=True)


if __name__ == '__main__':

    # enss = ['Ens{}'.format(i) for i in np.arange(1, 51, 1)]
    enss = ['Ens101']
    context = 'CH'

    log = snmcseq_utils.create_logger()
    for ens in enss:
        try:
            upload_marker_genes(ens, context, database=DATABASE)
            log.info('{} done!'.format(ens))
        except:
            log.info('{} skipped!'.format(ens))

