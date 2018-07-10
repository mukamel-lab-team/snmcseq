#!/usr/bin/env python3

import subprocess as sp

from __init__ import *
from snmcseq_utils import create_logger
import CEMBA_call_DMR

log = create_logger()

ens = 'Ens10'
allc_paths = glob.glob(os.path.join(PATH_ENSEMBLES, ens, 'allc_merged', 'allc_multimodal_v1_*.tsv'))
output_prefix = os.path.join(PATH_ENSEMBLES, ens, 'dmr', 'cgdmr_multimodal_v1')

CEMBA_call_DMR.call_DMR_wrapper(allc_paths, output_prefix, nprocs=8)
