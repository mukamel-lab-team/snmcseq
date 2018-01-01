#!/bin/bash

for perplexity in 10 15 20 25 30 50 100 150 200 250 500 1000
do
	input=./tsne/human_combined_hv1_hv2_CH_100000_out_v1_${perplexity}.tsv
	output=./tsne/plain/tsne_$(basename $input | cut -f 1 -d '.').pdf

	./plot_tsne_plain.py -i $input -o $output --title tsne_2human'_perplexity='${perplexity}
done 
