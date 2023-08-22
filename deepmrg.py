import numpy as np
import json
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale
from sklearn.metrics import classification_report
from tensorflow.keras.models import load_model
from Bio import SeqIO
import sys

#seed = 42
#np.random.seed(seed)
index_of_metals = {
    'Aluminium (Al)': 0,
    'Antimony (Sb)': 1,
    'Arsenic (As)': 2,
    'Bismuth (Bi)': 3,
    'Cadmium (Cd)': 4,
    'Chromium (Cr)': 5,
    'Cobalt (Co)': 6,
    'Copper (Cu)': 7,
    'Gallium (Ga)': 8,
    'Gold (Au)': 9,
    'Iron (Fe)': 10,
    'Lead (Pb)': 11,
    'Magnesium (Mg)': 12,
    'Manganese (Mn)': 13,
    'Mercury (Hg)': 14,
    'Molybdenum (Mo)': 15,
    'Nickel (Ni)': 16,
    'Selenium (Se)': 17,
    'Silver (Ag)': 18,
    'Tellurium (Te)': 19,
    'Tungsten (W)': 20,
    'Vanadium (V)': 21,
    'Zinc (Zn)': 22
}

index_of_metals_reverse = {
    0:'Aluminium (Al)',
    1:'Antimony (Sb)',
    2:'Arsenic (As)',
    3:'Bismuth (Bi)',
    4:'Cadmium (Cd)',
    5:'Chromium (Cr)',
    6:'Cobalt (Co)',
    7:'Copper (Cu)',
    8:'Gallium (Ga)',
    9:'Gold (Au)',
    10:'Iron (Fe)',
    11:'Lead (Pb)',
    12:'Magnesium (Mg)',
    13:'Manganese (Mn)',
    14:'Mercury (Hg)',
    15:'Molybdenum (Mo)',
    16:'Nickel (Ni)',
    17:'Selenium (Se)',
    18:'Silver (Ag)',
    19:'Tellurium (Te)',
    20:'Tungsten (W)',
    21:'Vanadium (V)',
    22:'Zinc (Zn)'
}
    
def compute_bitscore_feature_mat_with_clusters_of_ref_exp_and_annotation(bitscore_distribution, seq_fasta_file, ref_exp_mrg_clusters_dict, DeepMRG_model, output_prefix):

    x = []
    bitscore_if_no_alignment_found = 5 # defalut bitscore is set to 5 if no alignment report is found
    sorted_ref_exp_mrg_classes = sorted(ref_exp_mrg_clusters_dict.keys())

    index_of_proteins_in_list_x = {} # direct non-MRG will have index -1 and other seqs will have valid index as they will be feb to DeepMRG model.
    index_of_list_x = 0

    for sequence in SeqIO.parse(seq_fasta_file, "fasta"):

        prot_description = str(sequence.description)
        this_prot_id = prot_description

        # the case of proteins having no DIAMOND alignments at all with the 485 exp. verified MRG. This protein is annotated as non-MRG and not fed to DeepMRG model.
        if this_prot_id not in bitscore_distribution:
            index_of_proteins_in_list_x[this_prot_id] = -1
        
        else:

            this_prot_x = [] # this refers to current protein's bitscore distribution against 485 ref exp MRGs
            this_prot_bitscore_distribution = bitscore_distribution[this_prot_id] # this_prot_bitscore_distribution is a dict where key = ref_exp_gid and value = corresponding bitscore against the ref_exp_gid
            
            for i in sorted_ref_exp_mrg_classes:

                this_class_ref_exp_clusters = ref_exp_mrg_clusters_dict[i] # this_class_ref_exp_clusters is a dict where key = "Cluster p" (p is an int starting from 0) and value = ['BAC0036', 'BAC0025'].  
                n_clusters_of_this_class = len(this_class_ref_exp_clusters)
                sum_bitscore_of_this_class = 0

                for j in range(n_clusters_of_this_class):

                    cluster_name = "Cluster " + str(j)
                    sum_bitscore_of_this_cluster = 0
                    n_mrg_in_this_cluster = len(this_class_ref_exp_clusters[cluster_name])

                    for ref_exp_mrg in this_class_ref_exp_clusters[cluster_name]:

                        if ref_exp_mrg in this_prot_bitscore_distribution:
                            sum_bitscore_of_this_cluster+=this_prot_bitscore_distribution[ref_exp_mrg]
                        else:
                            sum_bitscore_of_this_cluster+=bitscore_if_no_alignment_found
                    avg_bitscore_of_this_cluster = sum_bitscore_of_this_cluster / n_mrg_in_this_cluster
                    sum_bitscore_of_this_class+=avg_bitscore_of_this_cluster
                
                avg_bitscore_of_this_class = sum_bitscore_of_this_class / n_clusters_of_this_class
                this_prot_x.append(avg_bitscore_of_this_class)

            x.append(this_prot_x)
            index_of_proteins_in_list_x[this_prot_id] = index_of_list_x
            index_of_list_x+=1

    if len(x): # at least one protein will be fed to DeepMRG
        x = np.array(x)
        testing_x_normalized = minmax_scale(x, axis=1) # sample (row) wise normalization instead of column wise
        model = load_model(DeepMRG_model)
        #print(model.summary())
        y_pred = model.predict(testing_x_normalized)
        
        confidence_threshold = 0.7
        #y_pred[y_pred >= confidence_threshold] = 1
        #y_pred[y_pred < confidence_threshold] = 0

    else: # all proteins have been discarded in the first filtering step
        pass


    deepmrg_annotation_filename = output_prefix + "_DeepMRG_annotation.tsv"

    annot_file = open(deepmrg_annotation_filename, "a")

    for sequence in SeqIO.parse(seq_fasta_file, "fasta"):

        prot_description = str(sequence.description)
        this_prot_id = prot_description

        if index_of_proteins_in_list_x[this_prot_id] == -1:
            annot_file.write(this_prot_id)
            annot_file.write("\t")
            annot_file.write("non-MRG")
            annot_file.write("\n")
        else:
            idx_of_this_prot_in_list_x = index_of_proteins_in_list_x[this_prot_id]
            deepmrg_pred_for_this_prot = list(y_pred[idx_of_this_prot_in_list_x])

            deepmrg_annot_for_this_prot = ""
            for i in range(len(deepmrg_pred_for_this_prot)):
                if deepmrg_pred_for_this_prot[i] >= confidence_threshold: # prediction confidence
                    deepmrg_annot_for_this_prot+=index_of_metals_reverse[i]
                    deepmrg_annot_for_this_prot+="(" + str(round(deepmrg_pred_for_this_prot[i]*100, 2)) + " %)"
                    deepmrg_annot_for_this_prot+=","
            
            if len(deepmrg_annot_for_this_prot):
                deepmrg_annot_for_this_prot = deepmrg_annot_for_this_prot[:-1]

                annot_file.write(this_prot_id)
                annot_file.write("\t")
                annot_file.write(deepmrg_annot_for_this_prot)
                annot_file.write("\n")

            else:
                annot_file.write(this_prot_id)
                annot_file.write("\t")
                annot_file.write("non-MRG")
                annot_file.write("\n")

    annot_file.close()

    return
    
if __name__ == "__main__":
    
    #ref_seq_path = '/home/muhitemon/DeepMRG_v3/processed_MRG-DB/MRG_DB_metal_wise_exp.fasta'
    testing_seq = str(sys.argv[1])
    testing_bitscore_distribution_file = str(sys.argv[2])
    clusters_of_ref_exp_mrg_json_file = str(sys.argv[3])
    DeepMRG_model = str(sys.argv[4])
    output_prefix = str(sys.argv[5])

    with open(clusters_of_ref_exp_mrg_json_file) as json_file:
        ref_exp_mrg_cluster_wise = json.load(json_file)
    
    with open(testing_bitscore_distribution_file) as json_file:
        testing_bitscore_distribution = json.load(json_file)
        
    compute_bitscore_feature_mat_with_clusters_of_ref_exp_and_annotation(testing_bitscore_distribution, testing_seq, ref_exp_mrg_cluster_wise, DeepMRG_model, output_prefix)
