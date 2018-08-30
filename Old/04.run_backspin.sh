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

# # human mch
# for n in 1000 2000 3000 4000 5000 6000
# do
# input='./data/binc/binc_mCH_human_combined_100000_summary_nmcc_v3.tsv'
# output=/cndd/fangming/snmcseq_dev/data/cluster/test_clusters/MB_v1_MB_EA_MB_EB/backspin_${n}cells.tsv
# min_var=0.15
# ./backspin.py -i $input -o $output -v $min_var -nsub $n 
# done

# human mch
for min_var in 0.63 0.65 0.68 
do
input='./data/binc/binc_mCH_human_combined_100000_summary_nmcc_v3.tsv.MB_v1'
	output=/cndd/fangming/snmcseq_dev/data/cluster/test_clusters/MB_v1/test_backspin_${min_var}p.tsv
./backspin.py -i $input -o $output -v $min_var 
done

# # human_hv1_hv2
# ./backspin.py -i ./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_backspin.tsv 

# # human_MB_EB
# ./backspin.py -i ./preprocessed/human_MB_EB_CH_100000_preproc.tsv -o ./preprocessed/human_MB_EB_CH_100000_backspin.tsv 

# # human_v1
# ./backspin.py -i ./preprocessed/human_v1_CH_100000_preproc.tsv -o ./preprocessed/human_v1_CH_100000_backspin.tsv 
