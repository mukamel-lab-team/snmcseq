#!/usr/bin/env python3

import pandas as pd
import os



# internal
def get_dmr_bed(input_f, outdir):
    """input_f: methylpy output (rms_results_collapsed.tsv)
    outdir: output directory
    """

    df = pd.read_table(input_f, index_col=['#chr', 'start', 'end'], dtype={'#chr': object})
    df_hypo = (df.loc[(df['number_of_dms']>=3) & (~df['hypomethylated_samples'].isnull()), 'hypomethylated_samples']
               .apply(lambda x: x.split(',')))

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        print('Created directory {}'.format(outdir))
        
    clsts = [i[len('methylation_level_'):] for i in df.filter(regex='^methylation_level_*').columns]
    for clst in clsts:
        df_hypo_cluster = df_hypo[df_hypo.apply(lambda x: (clst in x))]
        output = os.path.join(outdir, '{}.bed'.format(clst))
        df_hypo_cluster.to_csv(output, sep='\t', header=False, index=True, na_rep='NA')
        print("Saved to {}".format(output))
        print(df_hypo_cluster.shape)

    return 

if __name__ == '__main__':
    # get dmr bed files
        input_fs = ['/cndd/Public_Datasets/human_snmcseq/Ensembles/Ens5/dmr/cgdmr_human_mouse_v2-181120_human_v2_rms_results_collapsed.tsv', 
                '/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens0/dmr/cgdmr_human_mouse_v2-181120_mouse_rms_results_collapsed.tsv',
        ]
        for input_f in input_fs:
            if input_f.endswith('_rms_rsults_collapsed.tsv'):
                outdir = input_f[:-len('_rms_results_collapsed.tsv')] 
            get_dmr_bed(input_f, outdir)
