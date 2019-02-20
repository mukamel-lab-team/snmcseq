#!/usr/bin/env python3
"""allc_paths (the results from merging allc) set the format for downstream table and file names
allc_merged_mCG_$clustertype_$cluster_$ens.tsv.gz

- mc track files
mc_single_...$clustertype_$cluster_$ens_$chr.tsv
- mysql table
mc_single_..._$clustertype_$cluster_$ens_$chr
- track name 
mc_single_..._$clustertype_$cluster_$ens
- track vis name 
mCG_$cluster

"""

from __init__ import *

from natsort import natsorted
import subprocess as sp

from snmcseq_utils import create_logger
# import CEMBA_update_mysql
import CEMBA_merge_allc
import CEMBA_call_DMR
import CEMBA_get_dmr_bed
import CEMBA_allc_dmr_to_mc_single
import CEMBA_setup_annoj_v2
# from sqlalchemy import text 
from CEMBA_annoj_templates import *


# def upload_to_mysql_from_file(engine, file, table):
#     """
#     """
#     # create a tmp_table 
#     sqls = [
#         text('DROP TABLE IF EXISTS {}'.format(table)),
#         text('CREATE TABLE {} LIKE templates.mc_single'.format(table)),
#         text('LOAD DATA LOCAL INFILE \'{}\' INTO TABLE {}'.format(file, table)),
#         ]
#     for sql in sqls:
#         engine.execute(sql)
#     logging.info('Done loading {} to {}'.format(file, table))

#     return

def upload_to_mysql_from_file(database, file, table, script='./CEMBA_load_mc_track.sh'):
    """
    """
    logging.info('Loading {} to {}.{}'.format(file, database, table))
    cmd = '{} {} {} {}'.format(script, database, file, table)
    # print(cmd)
    sp.run(cmd, shell=True)

    return

if __name__ == '__main__':

    log = create_logger()

    enss = [
        'Ens38', #3F
        'Ens43', #1B
        'Ens44', #2E
        'Ens47', #4A
        'Ens50', #4E
        'Ens56', #1A
       #  'Ens57', #1C
        'Ens62', #4F
        'Ens63', #3E
        'Ens69', #2B
        'Ens70', #2D
        'Ens73', #4C
        'Ens74', #ACB
     ]  # example


    nprocs = 4
    # choose a clustering type
    cluster_type = 'cluster_mCHmCG_lv_npc50_k30'
    MERGE_ALLC = True
    CALL_DMR = True
    GET_DMR_BED = True
    GEN_MC_TRACKS = True
    UPLOAD_TO_MYSQL = True
    GEN_ANNOJ_TRACKS = True

    for ens in enss: 
        # ens = 'Ens1'
        # merge allc
        if MERGE_ALLC:
            CEMBA_merge_allc.merge_allc_CEMBA(ens, context='CG', cluster_type=cluster_type, database=DATABASE, 
                nprocs=nprocs,
                chunksize=1000000, n_chunk1=50, n_chunk2=20)

        # call DMRs
        allc_paths = natsorted(glob.glob(os.path.join(PATH_ENSEMBLES, ens, 
                        'allc_merged/allc_merged_mCG_{}_*_{}.tsv.gz'.format(cluster_type, ens)))) #natsorted
        dmr_dir = os.path.join(PATH_ENSEMBLES, ens, 'dmr')
        dmr_prefix = os.path.join(dmr_dir, 'dmr_allc_merged_mCG_{}'.format(cluster_type))
        if CALL_DMR:
            if not os.path.isdir(dmr_dir):
                os.makedirs(dmr_dir)
                logging.info("Created path: {}".format(dmr_dir))
            CEMBA_call_DMR.call_DMR_wrapper(allc_paths, dmr_prefix, nprocs=nprocs)

        # generate DMR bed files for each cluster
        dmr_summary_file = "{}_rms_results_collapsed.tsv".format(dmr_prefix)
        dmr_bed_dir = dmr_prefix
        if GET_DMR_BED:
            CEMBA_get_dmr_bed.get_dmr_bed(dmr_summary_file, dmr_bed_dir)

        # upload data to mysql
        dmr_files = natsorted(glob.glob(os.path.join(dmr_bed_dir, '*.bed')))  #natsorted
        annoj_dir = os.path.join(PATH_ENSEMBLES, ens, 'annoj') 
        if GEN_MC_TRACKS:
            # mc dmr to mc_single format
            use_dmrs = True
            if not os.path.isdir(annoj_dir):
                os.makedirs(annoj_dir)
            CEMBA_allc_dmr_to_mc_single.allc_dmr_to_mc_single(allc_paths, dmr_files, annoj_dir, use_dmrs) # make sure allc and dmr files corresponds

        # upload to mysql
        mc_track_files = natsorted(glob.glob(os.path.join(annoj_dir, 'mc_single_*.tsv')))
        mc_track_tables = [file.split('/')[-1][:-len('.tsv')] for file in mc_track_files]
        if UPLOAD_TO_MYSQL:
            for file, table in zip(mc_track_files, mc_track_tables):
                upload_to_mysql_from_file(DATABASE_ANNOJ, file, table)

        # gen browser tracks
        # title and output
        title = '{}_{}'.format(ens, cluster_type)
        track_names = natsorted(
                np.unique(['_'.join(table.split('_')[:-1]) for table in mc_track_tables]) # mysql table table_name_#chr 
                )
        track_vis_names = ['mCG_{}'.format(track_name.split('_')[-2]) for track_name in track_names] # ..._$clustertype_$cluster_$ens 
        output_php_dir = os.path.join(PATH_ENSEMBLES, ens, 'annoj_brainome') 
        output_html = os.path.join(output_php_dir, 'index_{}.html'.format(ens))
        if GEN_ANNOJ_TRACKS:
            # setup annoj
            mc_single_track_info = TrackInfo(
                TEMPLATE_MC_SINGLE_PHP, 
                TEMPLATE_INDEX_MC_SINGLE, 
                track_names, 
                track_vis_names, 
                scales='default', 
                use_scale=False,
                track_comment='mc_single',
                )

            # all tracks
            track_info_list = [
                mc_single_track_info,
            ]

            # active track names
            n = 5
            active_track_list = np.concatenate((
                mc_single_track_info.track_names[:n],
            ))

            CEMBA_setup_annoj_v2.setup_annoj_worker(title, track_info_list, output_php_dir, output_html, 
                        active_track_list=active_track_list)

            # transfer data to brainome (need ssh-keygen from banjo to brainome)
            brainome_html_dir = 'f7xie@brainome:/var/www/html/annoj_private/CEMBA'
            brainome_php_dir = 'f7xie@brainome:/var/www/html/annoj_private/CEMBA/browser/fetchers/mc_single'
            # output_html -> brainome_html
            cmd = 'rsync -vh --ignore-existing {} {}'.format(output_html, brainome_html_dir)
            print(cmd)
            sp.run(cmd, shell=True)

            # output_php_dir/php_mc_single/* -> brainome_php_dir 
            php_files = glob.glob(os.path.join(output_php_dir, 'php_mc_single', '*'))
            for _file in php_files:
                cmd = 'rsync -vh --ignore-existing {} {}'.format(_file, brainome_php_dir)
                print(cmd)
                sp.run(cmd, shell=True)


        # break
