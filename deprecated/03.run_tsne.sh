#!/bin/bash

for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
do
./preproc_and_TSNE.py -i ./tsne/human_combined_hv1_hv2_CH_100000.tsv -o ./tsne/human_combined_hv1_hv2_CH_100000_out_v1_${perplexity}.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv --normalize --perplexity ${perplexity} 
done

for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
do
./preproc_and_TSNE.py -i ./tsne/human_combined_MB_EB_CH_100000.tsv -o ./tsne/human_combined_MB_EB_CH_100000_out_v1_${perplexity}.tsv -m ./metadata/MB_EB_metadata_cells.tsv --normalize --perplexity ${perplexity} 
done

for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
do
./preproc_and_TSNE.py -i ./tsne/human_combined_v1_CH_100000.tsv -o ./tsne/human_combined_v1_CH_100000_out_v1_${perplexity}.tsv -m ./metadata/human_v1_metadata_cells.tsv --normalize --perplexity ${perplexity} 
done
