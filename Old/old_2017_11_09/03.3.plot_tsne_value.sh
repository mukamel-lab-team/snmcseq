#!/bin/bash

# global mCG human 1 and 2
./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv -l 'mCG/CG' -o ./tsne/colors/tsne_human_hv1_hv2_global_mCG.pdf -t 'human MFG and BA10 global mCG' -s

# global mCH human 1 and 2
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv -l mCH/CH -o ./tsne/colors/tsne_human_hv1_hv2_global_mCH.pdf -t 'human_1and2 global_mCH/CH' -s
# global mCH human 2
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_MB_EB_CH_100000_out_v1_25.tsv -m ./metadata/MB_EB_metadata_cells.tsv -l mCH/CH -o ./tsne/colors/tsne_human_MB_EB_global_mCH.pdf -t 'human_2 (MB_EB) global_mCH/CH' -s
# global mCH human 1
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_v1_CH_100000_out_v1_25.tsv -m ./metadata/human_v1_metadata_cells.tsv -l mCH/CH -o ./tsne/colors/tsne_human_v1_global_mCH.pdf -t 'human_1 global_mCH/CH' -s

# Genome coverage human 1 and 2 
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv -l '% Genome covered' -o ./tsne/colors/tsne_human_hv1_hv2_Genome_coverage.pdf -t 'human_1and2 genome coverage (%)' -s
# Genome coverage human 1 
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_v1_CH_100000_out_v1_25.tsv -m ./metadata/human_v1_metadata_cells.tsv -l 'Coverage (%)' -o ./tsne/colors/tsne_human_v1_Genome_coverage.pdf -t 'human_1 genome coverage (%)' -s
# Genome coverage human 2 
# ./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_MB_EB_CH_100000_out_v1_25.tsv -m ./metadata/MB_EB_metadata_cells.tsv -l '% Genome covered' -o ./tsne/colors/tsne_human_MB_EB_Genome_coverage.pdf -t 'human_2 genome coverage (%)' -s
