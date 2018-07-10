#!/bin/bash

input='/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_3C_171206/allc/allc_171213_CEMBA_mm_P56_P63_3C_MOp_CEMBA171206_3C_4_CEMBA171206_3C_5_H9_AD008_indexed.tsv.bgz'
output='/cndd/fangming/CEMBA/snmcseq_dev/test_mc_region_level.tsv' 
# bed_file='/cndd/kkolodzi/mouse_DMR/mouse_all_clusters_dmr_hyp.bed'
bed_file='/cndd/Public_Datasets/single_cell_methylome/DMRs/20170214/mouse/mouse_all_cluster_dms2_hypo.bed'
contexts='CH CG'

./CEMBA_mc_region_level.py -i $input -o $output -b $bed_file -c $contexts