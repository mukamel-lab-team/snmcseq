#!/bin/bash

# mCH
# input='./data/binc/binc_mCH_human_combined_100000_summary_nmcc_v3.tsv'
# output='./data/cluster/cluster_MB_v1_MB_EA_MB_EB/clusters_v1_binc_mCG_louvain.tsv' 
# ./louvain_jaccard.py -i $input -o $output

# mCG
input='./data/binc/binc_mCG_human_combined_100000_summary_nmcc_v3.tsv'
output='./data/cluster/cluster_MB_v1_MB_EA_MB_EB/clusters_v2_binc_mCG_louvain.tsv' 
./louvain_jaccard.py -i $input -o $output


# mCH and mCG
input='./data/binc/binc_mCHmCG_human_combined_100000_summary_nmcc_v3.tsv'
output='./data/cluster/cluster_MB_v1_MB_EA_MB_EB/clusters_v3_binc_mCHmCG_louvain.tsv' 
./louvain_jaccard.py -i $input -o $output



# # import tSNE coords
# tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCH_human_combined_100000_summary_nmcc_v3.tsv'
# tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCG_human_combined_100000_summary_nmcc_v3.tsv'
# tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCHmCG_human_combined_100000_summary_nmcc_v3.tsv'


