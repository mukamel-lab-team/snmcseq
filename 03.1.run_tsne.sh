#!/bin/bash

# # MB_v1 MB_EA MB_EB  mCHmCG bins 
# for perplexity in 10 15 20 25 30 50 100 150 200 250 300 500
# do
# ./tsne.py -i ./data/binc/binc_mCHmCG_human_combined_100000_summary_nmcc_v2.tsv -o ./data/tsne/tsne_perp${perplexity}_binc_mCHmCG_human_combined_100000_summary_nmcc_v2.tsv --perplexity ${perplexity} 
# done

# 
# MB_v1 MB_EA MB_EB  mCG bins 
for perplexity in 10 15 20 25 30 50 100 150 200 250 300 500
do
./tsne.py -i ./data/binc/binc_mCG_human_combined_100000_summary_nmcc_v2.tsv -o ./data/tsne/tsne_perp${perplexity}_binc_mCG_human_combined_100000_summary_nmcc_v2.tsv --perplexity ${perplexity} 
done

# # MB_v1 MB_EA MB_EB  mCH bins 
# for perplexity in 10 15 20 25 30 50 100 150 200 250 300 500
# do
# ./tsne.py -i ./data/binc/binc_mCH_human_combined_100000_summary_nmcc_v2.tsv -o ./data/tsne/tsne_perp${perplexity}_binc_mCH_human_combined_100000_summary_nmcc_v2.tsv --perplexity ${perplexity} 
# done
# 

# hs1-hs2 mCG gene

# for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
# do
# ./tsne_2.py -i ./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_out_v1_${perplexity}.tsv --perplexity ${perplexity} 
# done
# 
# 
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
