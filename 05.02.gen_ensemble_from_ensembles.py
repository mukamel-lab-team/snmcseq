#!/usr/bin/env python3
import pandas as pd
import CEMBA_init_ensemble_v2 
import snmcseq_utils


if __name__ == '__main__':

    log = snmcseq_utils.create_logger()

    f = '/cndd/fangming/CEMBA/data/19q2_ens.txt'
    df = pd.read_csv(f, sep='\t', header=None).values
    datasets = df[:, 0]
    ensid_list = df[:, 1]


    ens_id = 150
    ens_name = 'CEMBA_19Q2'
    ens_description = 'all CEMBA samples collected so far on cndd portal by 04/02/2019'
    ensid_list = sorted(ensid_list)  
    print(ensid_list)

    CEMBA_init_ensemble_v2.init_ensemble_from_ensembles(ens_id, ens_name, ens_description, 
                                ensid_list)


    # f = '/cndd/fangming/CEMBA/data/isocortex_Mar14_inhi.txt'
    # df = pd.read_csv(f, sep='\t', header=None).values
    # cell_list = df[:, 0]
    # 

    # ens_id = 149
    # ens_name = 'CEMBA_isocortex_19q1_inhibitory'
    # ens_description = 'Inhibitory neurons from all isocortex CEMBA samples (defined by Allen Brain Institute) on cndd portal by 03/12/2019'
    # ensid_list = sorted(ensid_list)  
    # enscell_list = cell_list

    # print(enscell_list)

    # CEMBA_init_ensemble_v2.init_ensemble_from_ensembles(ens_id, ens_name, ens_description, 
    #                             ensid_list, enscell_list=enscell_list)
