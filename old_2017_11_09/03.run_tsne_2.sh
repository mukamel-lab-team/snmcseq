#!/bin/bash

# hs1-hs2 mCG gene

for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
do
./tsne_2.py -i ./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_out_v1_${perplexity}.tsv --perplexity ${perplexity} 
done


# hv1_hv2
# for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
# do
# ./tsne_2.py -i ./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_out_v1_${perplexity}.tsv --perplexity ${perplexity} 
# done

# MB_EB
# for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
# do
# ./tsne_2.py -i ./preprocessed/human_MB_EB_CH_100000_preproc.tsv -o ./preprocessed/human_MB_EB_CH_100000_out_v1_${perplexity}.tsv --perplexity ${perplexity} 
# done
# 
# # human_v1
# for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
# do
# ./tsne_2.py -i ./preprocessed/human_v1_CH_100000_preproc.tsv -o ./preprocessed/human_v1_CH_100000_out_v1_${perplexity}.tsv --perplexity ${perplexity} 
# done
