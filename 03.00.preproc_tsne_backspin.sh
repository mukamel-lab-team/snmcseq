#!/bin/bash

input=./tsne/human_combined_hv1_hv2_CH_100000.tsv
output=./preprocessed/human_hv1_hv2_CH_100000_preproc.tsv
metadata=./metadata/human_hv1_hv2_metadata_cells.tsv
./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize

input=./tsne/human_combined_MB_EB_CH_100000.tsv
output=./preprocessed/human_MB_EB_CH_100000_preproc.tsv
metadata=./metadata/MB_EB_metadata_cells.tsv
./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize

input=./tsne/human_combined_v1_CH_100000.tsv
output=./preprocessed/human_v1_CH_100000_preproc.tsv
metadata=./metadata/human_v1_metadata_cells.tsv
./preproc_tsne_and_backspin.py -i $input -o $output -m $metadata --normalize
