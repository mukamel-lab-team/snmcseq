#!/usr/bin/env python3
"""Setup php and html files for AnnoJ browser

# prepare files
setup_annoj_worker(title, track_info_list, output_php_dir, output_html, 
    active_track_list=active_track_list)

Expected track_name:
atac_multimodal_v2_C123
rna_multimodal_v1_C1
mc_single_multimodal_v1_C1
...

TrackInfo example:
    mc_single_track_info = TrackInfo(
        TEMPLATE_MC_SINGLE_PHP, 
        TEMPLATE_INDEX_MC_SINGLE, 
        track_names, 
        track_vis_names, 
        scales='default', 
        use_scale=False,
        track_comment='mc_single',
        )


Fangming Xie
"""
import sys
sys.path.insert(0, '/cndd/fangming/CEMBA/snmcseq_dev')
from __init__ import *
from natsort import natsorted

import snmcseq_utils
from CEMBA_annoj_templates import *


if __name__ == '__main__':

    log = snmcseq_utils.create_logger()

    # title and output
    title = 'Ens10_multimodal_v2'
    output_php_dir = '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/annoj'
    output_html = (
        '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/annoj/'
        'index_{}.html'.format(title)
        )

    # annot file
    f = '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/annot_multimodal_v2.tsv'
    df_annot = pd.read_table(f, index_col='cluster_ID').iloc[:, 0]

    # track info list and active track list
    # mc_single tracks
    track_names = natsorted(np.ravel(pd.read_table(
        '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/mc_tracks/04.track_names.txt', 
        header=None).values)) 
    track_vis_names = ['{}_{}'.format(
        track_name.replace('multimodal_v2_', '').replace('mc_single', 'mCG'), 
        df_annot.loc[track_name.split('_')[4].replace('C', 'cluster_')], # get annot 
        ) for track_name in track_names] 

    mc_single_track_info = TrackInfo(
        TEMPLATE_MC_SINGLE_PHP, 
        TEMPLATE_INDEX_MC_SINGLE, 
        track_names, 
        track_vis_names, 
        scales='default', 
        use_scale=False,
        track_comment='mc_single',
        )
    # print(mc_single_track_info)

    # atac tracks 
    track_info = pd.read_table(
        '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/atac_clsts/04.track_names.txt',
        header=None, sep=' ', index_col=0, usecols=[0, 1])
    track_info = track_info.loc[natsorted(track_info.index.values), 1].reset_index()

    track_names = track_info[0].values
    scales = track_info[1].values
    scales = scales/np.mean(scales)
    track_vis_names = ['{}_{}'.format(
        track_name.replace('multimodal_v2_', '').upper(), 
        df_annot.loc[track_name.split('_')[3].replace('C', 'cluster_')], # get annot 
        ) for track_name in track_names] 
    # print(track_vis_names)

    atac_track_info = TrackInfo(
        TEMPLATE_ATAC_PHP, 
        TEMPLATE_INDEX_ATAC, 
        track_names, 
        track_vis_names, 
        scales=scales, 
        use_scale=True, # has to be True
        track_comment='atac')

    # rna tracks 
    track_info = pd.read_table(
        '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/snrna_clsts/04.track_names.txt',
        header=None, sep=' ', index_col=0, usecols=[0, 1])
    track_info = track_info.loc[natsorted(track_info.index.values), 1].reset_index()

    track_names = track_info[0].values
    scales = track_info[1].values
    scales = scales/np.mean(scales)
    track_vis_names = ['{}_{}'.format(
        track_name.replace('multimodal_v2_', '').upper(), 
        df_annot.loc[track_name.split('_')[3].replace('C', 'cluster_')], # get annot 
        ) for track_name in track_names] 
    # print(track_vis_names)

    rna_track_info = TrackInfo(
        TEMPLATE_RNA_PHP, 
        TEMPLATE_INDEX_RNA, 
        track_names, 
        track_vis_names, 
        scales=scales, 
        use_scale=True, # has to be True 
        track_comment='rna')

    # enhancer tracks
    track_names = natsorted(np.ravel(pd.read_table(
        '/cndd/fangming/CEMBA/annoj_browser/multimodal_v2/reptile_enhancer_tracks/04.track_names.txt',
        header=None).values)) 
    track_vis_names = track_names 

    enhancer_track_info = TrackInfo(
        TEMPLATE_BED_PHP, 
        TEMPLATE_INDEX_DMR, 
        track_names, 
        track_vis_names, 
        scales='default', 
        use_scale=False, 
        track_comment='enhancer')

    # all tracks
    track_info_list = [
        mc_single_track_info,
        atac_track_info,
        rna_track_info,
        enhancer_track_info,
    ]

    # active_track_list = []
    n = 5
    active_track_list = np.concatenate((
        mc_single_track_info.track_names[:n],
        atac_track_info.track_names[:n],
        rna_track_info.track_names[:n],
        enhancer_track_info.track_names,
    ))

    # prepare files
    setup_annoj_worker(title, track_info_list, output_php_dir, output_html, 
        active_track_list=active_track_list)

