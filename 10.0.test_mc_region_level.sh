#!/bin/bash

input='/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens10/allc_merged/allc_multimodal_v2_clst10.tsv.gz'
output='/cndd/fangming/CEMBA/snmcseq_dev/test_mc_region_level.tsv' 
# bed_file='/cndd/kkolodzi/mouse_DMR/mouse_all_clusters_dmr_hyp.bed'
bed_file='/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens100/dmr/cgdmr_multimodal_v2_rms_results_collapsed.tsv.DMR_3dms.bed'
# contexts='CH CG'
contexts='CG'

./CEMBA_mc_region_level.py -i $input -o $output -b $bed_file -c $contexts