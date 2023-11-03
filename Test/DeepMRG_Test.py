import numpy as np
import json
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale
from sklearn.metrics import classification_report
from tensorflow.keras.models import load_model
from Bio import SeqIO
import joblib
import argparse
import os
import sys

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

def create_lablels_of_protein_seq(metals):
    '''
    argument metals can be like Zinc (Zn), Copper (Cu) i.e., csv
    '''
    label = [0] * 23 # there are 23 types of metals in BacMet DB
    metal_list = metals.split(',')
    for i in metal_list:
        if i.strip() in index_of_metals:
            label[index_of_metals[i.strip()]] = 1
    
    return label
    
def compute_bitscore_feature_mat_with_clusters_of_ref_exp_and_label(bitscore_distribution, seq_fasta_file, ref_exp_mrg_clusters_dict):

    x = []
    y = []
    bitscore_if_no_alignment_found = 0 # defalut bitscore is set to 0 if no alignment report is found
    sorted_ref_exp_mrg_classes = sorted(ref_exp_mrg_clusters_dict.keys())

    for sequence in SeqIO.parse(seq_fasta_file, "fasta"):

        prot_description = str(sequence.description)
        this_prot_id = prot_description

        this_prot_x = [] # this refers to current protein's bitscore distribution against 485 ref exp MRGs

        # the case of proteins having no DIAMOND alignments at all with the 485 exp. verified MRG
        if this_prot_id not in bitscore_distribution:

            for i in sorted_ref_exp_mrg_classes:
                
                this_prot_x.append(bitscore_if_no_alignment_found)
        
        else:

            this_prot_bitscore_distribution = bitscore_distribution[this_prot_id] # this_prot_bitscore_distribution is a dict where key = ref_exp_gid and value = corresponding bitscore against the ref_exp_gid
            
            for i in sorted_ref_exp_mrg_classes:
    
                this_class_ref_exp_clusters = ref_exp_mrg_clusters_dict[i] # this_class_ref_exp_clusters is a dict where key = "Cluster p" (p is an int starting from 0) and value = ['BAC0036', 'BAC0025'].  
                n_clusters_of_this_class = len(this_class_ref_exp_clusters)
                sum_bitscore_of_this_class = 0
                max_bitscore_of_this_class = 0

                for j in range(n_clusters_of_this_class):

                    cluster_name = "Cluster " + str(j)
                    sum_bitscore_of_this_cluster = 0
                    max_bitscore_of_this_cluster = 0
                    n_mrg_in_this_cluster = len(this_class_ref_exp_clusters[cluster_name])

                    for ref_exp_mrg in this_class_ref_exp_clusters[cluster_name]:

                        if ref_exp_mrg in this_prot_bitscore_distribution:
                            #sum_bitscore_of_this_cluster+=this_prot_bitscore_distribution[ref_exp_mrg]
                            if this_prot_bitscore_distribution[ref_exp_mrg] > max_bitscore_of_this_cluster:
                                max_bitscore_of_this_cluster = this_prot_bitscore_distribution[ref_exp_mrg]
                        else:
                            #sum_bitscore_of_this_cluster+=bitscore_if_no_alignment_found
                            if bitscore_if_no_alignment_found > max_bitscore_of_this_cluster:
                                max_bitscore_of_this_cluster = bitscore_if_no_alignment_found

                    #avg_bitscore_of_this_cluster = sum_bitscore_of_this_cluster / n_mrg_in_this_cluster
                    #sum_bitscore_of_this_class+=avg_bitscore_of_this_cluster
                    #sum_bitscore_of_this_class+=sum_bitscore_of_this_cluster
                    max_bitscore_of_this_class = max(max_bitscore_of_this_class, max_bitscore_of_this_cluster)

                #avg_bitscore_of_this_class = sum_bitscore_of_this_class / n_clusters_of_this_class
                this_prot_x.append(max_bitscore_of_this_class)


        this_prot_label = prot_description.split('|')[-1]
        if this_prot_label.lower() == "multimetal":
            this_prot_label = prot_description.split('|')[-2]

        x.append(this_prot_x)
        y.append(create_lablels_of_protein_seq(this_prot_label))


    x = np.array(x)
    y = np.array(y)

    return x, y
    
if __name__ == "__main__":
    
    testing_set = str(sys.argv[1]) # Just a string. One of ["TEST", "LOW", "STM6070", "fold1", "fold2", "fold3", "fold4", "fold5"]
    testing_seq_path = str(sys.argv[2])
    testing_bitscore_distribution_file = str(sys.argv[3])
    clusters_of_ref_exp_mrg_json_file = str(sys.argv[4])
    DeepMRG_model1 = str(sys.argv[5])
    DeepMRG_model2 = str(sys.argv[6])
    DeepMRG_model3 = str(sys.argv[7])
    DeepMRG_model4 = str(sys.argv[8])
    DeepMRG_model5 = str(sys.argv[9])
    
    
    with open(clusters_of_ref_exp_mrg_json_file) as json_file:
        ref_exp_mrg_cluster_wise = json.load(json_file)
        
    
    with open(testing_bitscore_distribution_file) as json_file:
        testing_bitscore_distribution = json.load(json_file)
        
    testing_x, testing_y_true = compute_bitscore_feature_mat_with_clusters_of_ref_exp_and_label(testing_bitscore_distribution, testing_seq_path, ref_exp_mrg_clusters_dict=ref_exp_mrg_cluster_wise)

    
    testing_x_normalized = minmax_scale(testing_x, axis=1) # sample (row) wise normalization instead of column wise
    
    model1 = load_model(DeepMRG_model1)
    model2 = load_model(DeepMRG_model2)
    model3 = load_model(DeepMRG_model3)
    model4 = load_model(DeepMRG_model4)
    model5 = load_model(DeepMRG_model5)
    
    if testing_set == "TEST" or testing_set == "LOW" or testing_set == "STM6070":
        y_pred1 = model1.predict(testing_x_normalized)
        y_pred2 = model2.predict(testing_x_normalized)
        y_pred3 = model3.predict(testing_x_normalized)
        y_pred4 = model4.predict(testing_x_normalized)
        y_pred5 = model5.predict(testing_x_normalized)
        
        y_preds = [y_pred1, y_pred2, y_pred3, y_pred4, y_pred5]
        Y_PRED = np.sum(y_preds, axis=0)
        
        confidence_threshold = 3.5
        Y_PRED[Y_PRED < confidence_threshold] = 0
        Y_PRED[Y_PRED >= confidence_threshold] = 1
        
        report = classification_report(testing_y_true, Y_PRED, output_dict=False,
		        target_names=['Aluminium (Al)',
		                        'Antimony (Sb)',
		                        'Arsenic (As)',
		                        'Bismuth (Bi)',
		                        'Cadmium (Cd)',
		                        'Chromium (Cr)',
		                        'Cobalt (Co)',
		                        'Copper (Cu)',
		                        'Gallium (Ga)',
		                        'Gold (Au)',
		                        'Iron (Fe)',
		                        'Lead (Pb)',
		                        'Magnesium (Mg)',
		                        'Manganese (Mn)',
		                        'Mercury (Hg)',
		                        'Molybdenum (Mo)',
		                        'Nickel (Ni)',
		                        'Selenium (Se)',
		                        'Silver (Ag)',
		                        'Tellurium (Te)',
		                        'Tungsten (W)',
		                        'Vanadium (V)',
		                        'Zinc (Zn)'])
	
	# 5-fold cross-validation
    else:
        if testing_set == "fold1":
            y_pred = model1.predict(testing_x_normalized)
        elif testing_set == "fold2":
            y_pred = model2.predict(testing_x_normalized)
        elif testing_set == "fold3":
            y_pred = model3.predict(testing_x_normalized)
        elif testing_set == "fold4":
            y_pred = model4.predict(testing_x_normalized)
        elif testing_set == "fold5":
            y_pred = model5.predict(testing_x_normalized)
        
        confidence_threshold = 0.7
        y_pred[y_pred < confidence_threshold] = 0
        y_pred[y_pred >= confidence_threshold] = 1
        
        report = classification_report(testing_y_true, y_pred, output_dict=False,
		        target_names=['Aluminium (Al)',
		                        'Antimony (Sb)',
		                        'Arsenic (As)',
		                        'Bismuth (Bi)',
		                        'Cadmium (Cd)',
		                        'Chromium (Cr)',
		                        'Cobalt (Co)',
		                        'Copper (Cu)',
		                        'Gallium (Ga)',
		                        'Gold (Au)',
		                        'Iron (Fe)',
		                        'Lead (Pb)',
		                        'Magnesium (Mg)',
		                        'Manganese (Mn)',
		                        'Mercury (Hg)',
		                        'Molybdenum (Mo)',
		                        'Nickel (Ni)',
		                        'Selenium (Se)',
		                        'Silver (Ag)',
		                        'Tellurium (Te)',
		                        'Tungsten (W)',
		                        'Vanadium (V)',
		                        'Zinc (Zn)'])
		
    	                       
    output_file = testing_set + "_classification_report.txt"
    with open(output_file, "w") as file:
    	file.write(report)