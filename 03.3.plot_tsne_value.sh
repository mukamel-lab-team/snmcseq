# global mCH
./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_MB_EB_CH_100000_out_v1_25.tsv -m ./metadata/MB_EB_metadata_cells.tsv -l mCH/CH -o ./tsne/colors/tsne_human_MB_EB_global_mCH.pdf -t human_MB_EB_global_mCH/CH
# global mCH
./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_hv1_hv2_CH_100000_out_v1_25.tsv -m ./metadata/human_hv1_hv2_metadata_cells.tsv -l mCH/CH -o ./tsne/colors/tsne_human_hv1_hv2_global_mCH.pdf -t human_hv1_hv2_global_mCH/CH
# Genome coverage 
./scripts_plot/plot_tsne_value.py -i ./preprocessed/human_MB_EB_CH_100000_out_v1_25.tsv -m ./metadata/MB_EB_metadata_cells.tsv -l '% Genome covered' -o ./tsne/colors/tsne_human_MB_EB_Genome_coverage.pdf -t human_MB_EB_%_genome_coverage
