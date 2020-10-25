#!/bin/bash

input='/cndd2/fangming/scf_enhancers/data/bulk/round3/mc/Round3.1-1-1.CGN-Merge.allc.tsv.gz'
output='test.tsv' 
bed_file='/cndd2/fangming/scf_enhancers/scripts/enhancers/test.bed'
# contexts='CH CG'
contexts='CG'

./CEMBA_mc_region_level.py -i $input -o $output -b $bed_file -c $contexts
