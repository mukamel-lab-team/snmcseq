# Library with utilities from Chris Keown's "mypy"

def get_human_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,23)]
    if include_x:
        chromosomes.append('X')
    return chromosomes


def read_gencode_human(version='v19', pc=False):
	# pc = protein coding
    prefix = '/cndd/projects/Public_Datasets/references/hg19/transcriptome/'
    if pc:
        fname= prefix+'gencode.'+version+'.annotation_genes_pc_mypy.tsv'
    else:
        fname= prefix+'gencode.'+version+'.annotation_genes_mypy.tsv'
    return pd.read_csv(fname, sep="\t")

def tabix_summary(records, context="CH", cap=0):

    mc = 0
    c = 0

    if context == "CH":
        contexts = get_mCH_contexts()
    elif context == "CA":
        contexts = ["CAA","CAC","CAG","CAT"]
    elif context == "CT":
        contexts = ["CTA","CTC","CTG","CTT"]
    elif context == "CAG":
        contexts = ["CAG"]
    elif context == "CAC":
        contexts = ["CAC"]
    else:
        contexts = get_mCG_contexts()+['CGN']

    if cap > 0:
        for record in records:
            if record[3] in contexts:
                if int(record[5]) <= cap:
                    mc += int(record[4])
                    c += int(record[5])
    else:

        for record in records:
            if record[3] in contexts:
                mc += int(record[4])
                c += int(record[5])

    return mc, c


def get_mCH_contexts():
    contexts = []
    for base1 in ['A','C','T']:
        for base2 in ['A','C','G','T']:
            contexts.append('C' + base1 + base2)
    return contexts
