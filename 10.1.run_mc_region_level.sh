#!/bin/bash
# the bash script to call the following python script

# usage: CEMBA_run_mc_region_level.py [-h] -d DATASETS [DATASETS ...] -o
#                                     OUTPUT_DIRNAME -b BED_FILE
#                                     [-c CONTEXTS [CONTEXTS ...]] [-n NPROCS]

# optional arguments:
#   -h, --help            show this help message and exit
#   -d DATASETS [DATASETS ...], --datasets DATASETS [DATASETS ...]
#                         list of datasets
#   -o OUTPUT_DIRNAME, --output_dirname OUTPUT_DIRNAME
#                         output dirname UNDER dataset folders (example:
#                         gene_level, binc, atac_peaks, ...)
#   -b BED_FILE, --bed_file BED_FILE
#                         bed file
#   -c CONTEXTS [CONTEXTS ...], --contexts CONTEXTS [CONTEXTS ...]
#                         list of contexts: CH/CG/... default: CH CG CA
#   -n NPROCS, --nprocs NPROCS
#                         number of processes (default 1)

datasets="CEMBA_3C_171206 \
CEMBA_3C_171207 \
CEMBA_4B_171212 \
CEMBA_4B_171213 \
CEMBA_4B_180104"

outdirname="1kb_bins"
bed_file="/cndd/fangming/iGenome/mm10/mm10_1kb.bed"
contexts="CH CG"
nprocs=8

./CEMBA_run_mc_region_level.py -d $datasets -o $outdirname -b $bed_file -c $contexts -n $nprocs

