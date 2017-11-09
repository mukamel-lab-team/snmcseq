#!/usr/bin/env python3

import pandas as pd

input_fname = './metadata/human_v1_metadata_cells.tsv'

df = pd.read_table(input_fname, header=0)
df['Sample'] = [item+'_R1' for item in df['Sample'].tolist()]

df.to_csv(input_fname+'vvv', sep='\t', header=True, index=False)