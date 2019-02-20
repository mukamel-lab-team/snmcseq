#!/usr/bin/env python3

import subprocess as sp

from __init__ import *
from snmcseq_utils import create_logger
import CEMBA_call_DMR

log = create_logger()

ens = 'Ens5'
PATH_ENSEMBLES = '/cndd/Public_Datasets/human_snmcseq/Ensembles'
allc_paths = glob.glob(os.path.join(PATH_ENSEMBLES, ens, 'allc_merged', 'allc_human_mouse_v2-181120_*_human_v2.tsv.gz'))
output_prefix = os.path.join(PATH_ENSEMBLES, ens, 'dmr', 'cgdmr_human_mouse_v2-181120_human_v2')
file_suffix = '.tsv.gz'

CEMBA_call_DMR.call_DMR_wrapper(allc_paths, output_prefix, file_suffix=file_suffix, nprocs=8)
