import pandas as pd 
import json

def compute_bitscore(alignment):
    '''
    parameters:
        alignment is a panda dataframe
    returns:
        bitscore_distribution: dict where key is gene_id and value is a dict of ref_exp gene_id (key) and corresponding bitscore with this ref_exp gene_id (value)
    '''
    bitscore_distribution = {}
    for i in range(len(alignment)):

        qtitle = alignment.iloc[i]['q_mrg']
        split_qtitle = qtitle.split('|')
        qid = str(split_qtitle[0].strip())

        ref_title = alignment.iloc[i]['ref_mrg']
        split_ref_title = ref_title.split('|')
        ref_id = str(split_ref_title[0].strip())

        if qid not in bitscore_distribution:
            bitscore_distribution[qid] = {ref_id:float(alignment.iloc[i]['bitscore'])}
        else:
            bitscore_distribution[qid][ref_id] = float(alignment.iloc[i]['bitscore'])

    return bitscore_distribution

if __name__ == '__main__':

    alignment_file = '/home/muhitemon/DeepMRG/bitscore_distribution_with_ref_exp/test_diamond_with_ref_exp_from_single_linkage40.tsv'
    diamond_output = pd.read_csv(alignment_file, sep='\t', names=['q_mrg','ref_mrg', 'pident', 'bitscore', 'evalue', 'q_mrg_len', 'ref_mrg_len', 'q_mrg_start', 'q_mrg_end', 'ref_mrg_start', 'ref_mrg_end'], header=None)
    
    bitscore_distribution = compute_bitscore(diamond_output)

    print(len(bitscore_distribution))
    
    with open('test_bitscore_distribution_with_ref_exp_from_single_linkage40.json', 'w') as fp:
        json.dump(bitscore_distribution, fp)
    

    #'WP_001114232.1'