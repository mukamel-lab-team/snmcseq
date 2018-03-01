#!/bin/bash

# human mch and mcg
# input='./data/binc/binc_mCHmCG_human_combined_100000_summary_nmcc_v2.tsv'
# output='./data/backspin/clusters_binc_mCHmCG_human_combined_100000_summary_nmcc_v2.tsv'
# ./backspin.py -i $input -o $output 
# 
# # human mcg
# input='./data/binc/binc_mCG_human_combined_100000_summary_nmcc_v2.tsv'
# output='./data/backspin/clusters_binc_mCG_human_combined_100000_summary_nmcc_v2.tsv'
# ./backspin.py -i $input -o $output 
# 
# human mch
input='./data/binc/binc_mCH_human_combined_100000_summary_nmcc_v2.tsv'
output='./data/backspin/clusters_binc_mCH_human_combined_100000_summary_nmcc_v2.tsv'
./backspin.py -i $input -o $output 

# # human_hv1_hv2
# ./backspin.py -i ./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_backspin.tsv 

# # human_MB_EB
# ./backspin.py -i ./preprocessed/human_MB_EB_CH_100000_preproc.tsv -o ./preprocessed/human_MB_EB_CH_100000_backspin.tsv 

# # human_v1
# ./backspin.py -i ./preprocessed/human_v1_CH_100000_preproc.tsv -o ./preprocessed/human_v1_CH_100000_backspin.tsv 
