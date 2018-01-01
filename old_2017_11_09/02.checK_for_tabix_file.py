import mypy
import os
import multiprocessing as mp
import numpy as np



samples = mypy.get_single_cell_samples_human_FNL()

allc_dir = '/cndd/projects/Public_Datasets/single_cell_methylome/hs_20161229/'

# samples = samples[:2]

for sample in samples:
    for chromosome in mypy.get_human_chromosomes() + ['Y']:
        if not os.path.isfile(allc_dir+sample+'_bismark/allc_'+sample+'_'+chromosome+'.tsv.gz'):
            print('Doesnt exist: ' + allc_dir+sample+'_bismark/allc_'+sample+'_'+chromosome+'.tsv.gz')

# allc_'+sample+'_'+chromosome+'.tsv.gz')
#     if os.path.isdir(allc_dir+sample+'_bismark'):
#         print(allc_dir+sample+'_bismark')
#         os.chdir(allc_dir+sample+'_bismark')
#         for chromosome in mypy.get_human_chromosomes() + ['Y']:
#             # os.system('bgzip -c allc_'+sample+'_'+chromosome+'.tsv > allc_'+sample+'_'+chromosome+'.tsv.gz')
#             # print('allc_'+sample+'_'+chromosome+'.tsv.gz')
#             os.system('tabix -f -s 1 -b 2 -e 2 -S 1 allc_'+sample+'_'+chromosome+'.tsv.gz')




# def process(samples):
#     for sample in samples:
#         if os.path.isdir(allc_dir+sample+'_bismark'):
#             print(allc_dir+sample+'_bismark')
#             os.chdir(allc_dir+sample+'_bismark')
#             for chromosome in mypy.get_human_chromosomes() + ['Y']:
#                 # os.system('bgzip -c allc_'+sample+'_'+chromosome+'.tsv > allc_'+sample+'_'+chromosome+'.tsv.gz')
#                 # print('allc_'+sample+'_'+chromosome+'.tsv.gz')
#                 os.system('tabix -f -s 1 -b 2 -e 2 -S 1 allc_'+sample+'_'+chromosome+'.tsv.gz')
#         else:
#             print('WARNING: THE BELOW DIRECTORY DOES NOT EXIST.')
#             print(allc_dir+sample+'_bismark')

#     return True


# procs = min(len(samples), 16)
# p = mp.Pool(processes=procs)
# split_samples = np.array_split(samples,procs)
# pool_results = p.map(process, split_samples)
# p.close()
# p.join()

