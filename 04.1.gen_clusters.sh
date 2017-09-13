#!/bin/bash
# hv1_hv2 
./postproc_backspin.py -i ./preprocessed/human_hv1_hv2_CH_100000_backspin.tsv -o ./preprocessed/human_hv1_hv2_CH_100000_clusters.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv

# hv1 
./postproc_backspin.py -i ./preprocessed/human_v1_CH_100000_backspin.tsv -o ./preprocessed/human_v1_CH_100000_clusters.tsv -m ./metadata/human_v1_metadata_cells.tsv

# MB_EB
./postproc_backspin.py -i ./preprocessed/human_MB_EB_CH_100000_backspin.tsv -o ./preprocessed/human_MB_EB_CH_100000_clusters.tsv -m ./metadata/MB_EB_metadata_cells.tsv

