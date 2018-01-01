#!/bin/bash

# joint clustering 
# ./scripts_plot/plot_tsne_2labels.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv -r ./preprocessed/human_hv1_hv2_CH_100000_clusters.tsv -l 'cluster_ID' --merge_col 'sample' -o ./preprocessed/human_hv1_hv2_CH_100000_tsne_clusters.pdf --title 'human 1 and 2 joint backspin clustering (82 clusters)' --legend_mode -1 -s

# joint clustering 3d
# ./scripts_plot/plot_tsne_2labels_3d.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_perp25_ntSNE3.tsv -r ./preprocessed/human_hv1_hv2_CH_100000_clusters.tsv -l 'cluster_ID' --merge_col 'sample' -o ./preprocessed/human_hv1_hv2_CH_100000_tsne_clusters.pdf --title 'human 1 and 2 joint backspin clustering (82 clusters)' --legend_mode -1 -s

# human2 clustering 
# ./plot_tsne_clusters.py -i ./preprocessed/human_MB_EB_CH_100000_out_v1_25.tsv -r ./preprocessed/human_MB_EB_CH_100000_clusters.tsv -l 'cluster_ID' --merge_col 'sample' -o ./preprocessed/human_MB_EB_CH_100000_tsne_clusters.pdf -t 'human MB EB (human 2) backspin clusters (42 clusters)' --legend_mode 1 -s

# human1 clustering 
# ./plot_tsne_clusters.py -i ./preprocessed/human_v1_CH_100000_out_v1_25.tsv -r ./preprocessed/human_v1_CH_100000_clusters.tsv -l 'cluster_ID' --merge_col 'sample' -o ./preprocessed/human_v1_CH_100000_tsne_clusters.pdf -t 'human 1 backspin clusters (50 clusters)' --legend_mode 1 -s

# joint clustering label clusters
./scripts_plot/plot_tsne_2labels.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv -r ./metadata/human_hv1_hv2_metadata_cells.tsv -l 'cluster_label' --merge_col 'Sample' -o ./preprocessed/human_hv1_hv2_CH_100000_tsne_cell_types.pdf --title 'human 1 and 2 tsne labeled with cell types' -s
