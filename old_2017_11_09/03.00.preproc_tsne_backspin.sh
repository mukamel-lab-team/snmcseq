#!/bin/bash


# hs1-hs2-mch genebody
input=./metadata/genebody_mCH_combined-hs1-hs2.mat
num_ex_cols=1 # first 1 column is gene id
context=CH
output=./preprocessed/human_hv1_hv2_mCH_genebody_nmcc.mat
metadata=./metadata/human_hv1_hv2_metadata_cells.tsv
base_call_cutoff=100
./preproc_tsne_and_backspin.py -i $input -ex $num_ex_cols -cx $context -o $output -m $metadata --normalize -b $base_call_cutoff


# # hs1-hs2-mcg genebody
# input=./metadata/genebody_mCG_combined-hs1-hs2.mat
# num_ex_cols=1 # first 1 column is gene id
# context=CG
# output=./preprocessed/human_hv1_hv2_CG_genebody_nmcc.mat
# metadata=./metadata/human_hv1_hv2_metadata_cells.tsv
# base_call_cutoff=25
# ./preproc_tsne_and_backspin.py -i $input -ex $num_ex_cols -cx $context -o $output -m $metadata --normalize -b $base_call_cutoff
# 

# input=./tsne/human_combined_hv1_hv2_CH_100000.tsv
# output=./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv
# metadata=./metadata/human_hv1_hv2_metadata_cells.tsv
# ./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize
# 
# input=./tsne/human_combined_MB_EB_CH_100000.tsv
# output=./preprocessed/human_MB_EB_CH_100000_preproc.tsv
# metadata=./metadata/MB_EB_metadata_cells.tsv
# ./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize
# 
# input=./tsne/human_combined_v1_CH_100000.tsv
# output=./preprocessed/human_v1_CH_100000_preproc.tsv
# metadata=./metadata/human_v1_metadata_cells.tsv
# ./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize
