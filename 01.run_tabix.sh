#!/bin/bash
#
# Loop through all allc files and run tabix indexing

allc_dir=$1

find $allc_dir -type f -name 'allc*.tsv.gz' | xargs -P16 -I{} allc_tabix {}
