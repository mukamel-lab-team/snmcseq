#!/bin/bash

# mCG binc 100kb 
input=./data/binc/binc_mCG_human_combined_100000_summary.tsv
num_ex_cols=2 # first 2 columns are chr and bin 
context=CG
output=./data/binc/binc_mCG_human_combined_100000_summary_normalized.tsv
metadata=./data/metadata/metadata_human_combined_updated.tsv
base_call_cutoff=30
sufficient_coverage_fraction=0.98
./preproc_tsne_and_backspin.py -i $input -ex $num_ex_cols -cx $context -o $output -m $metadata --normalize -b $base_call_cutoff -s $sufficient_coverage_fraction

# # mCH binc 100kb 
# input=./data/binc/binc_mCH_human_combined_100000_summary.tsv
# num_ex_cols=2 # first 2 columns are chr and bin 
# context=CH
# output=./data/binc/binc_mCH_human_combined_100000_summary_normalized.tsv
# metadata=./data/metadata/metadata_human_combined_updated.tsv
# base_call_cutoff=100
# ./preproc_tsne_and_backspin.py -i $input -ex $num_ex_cols -cx $context -o $output -m $metadata --normalize -b $base_call_cutoff
# 
