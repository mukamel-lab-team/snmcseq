#!/bin/bash

input='/cndd/Public_Datasets/single_cell_methylome/mapped_mouse_161224_q30_allc/Pool_1168_AD002_indexed'
output='./test_mc_region_level_old.tsv'
bed_file='/cndd/junhao/genomes/hg19/genomicRegions/UCSC_repeatMasker_LINE.bed'
contexts='CH CG'

./mc_region_level.py -i $input -o $output -b $bed_file -c $contexts