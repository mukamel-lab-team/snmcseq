#!/bin/bash

./scripts_plot/plot_tsne_2labels.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv -l 'Library pool' -s -r ./preprocessed/human_hv1_hv2_CH_100000_clusters.tsv -l2 cluster_ID --merge_col2 sample
