import h5py
import numpy as np
import json
import tensorflow as tf
from sklearn.utils import shuffle
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Input
from tensorflow.keras.models import Sequential, load_model, save_model, Model
from tensorflow.keras.layers import BatchNormalization, Dropout
from Bio import SeqIO

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

def create_lablels_of_protein_seq(metals):
    '''
    argument metals can be like Zinc (Zn), Copper (Cu) i.e., csv
    '''
    label = [0] * 23 # there are 23 types of metals in BacMet DB
    metal_list = metals.split(',')
    for i in metal_list:
        label[index_of_metals[i.strip()]] = 1
    
    return label

def find_ref_exp_mrg_class(metals):
    '''
    argument metals can be like Zinc (Zn), Copper (Cu) i.e., csv
    returns the metal class in sorted way like Antimony (Sb),Arsenic (As)
    '''
    L = []
    ref_exp_mrg_class = ""
    metal_list = metals.split(',')
    for i in metal_list:
        L.append((index_of_metals[i.strip()], i.strip()))
    
    L.sort()
    for i, metal in L:
        ref_exp_mrg_class+=metal + ','
    ref_exp_mrg_class = ref_exp_mrg_class[:-1] # omit the last comma
    return ref_exp_mrg_class

def get_model(number_of_neurons_in_each_hlayer = [50, 40, 30]):

    dimension_of_feature = 66

    input1 = Input(shape=(dimension_of_feature,))
    
    
    for i in range(len(number_of_neurons_in_each_hlayer)):
        if i == 0:
            x = Dense(number_of_neurons_in_each_hlayer[i], activation= 'elu')(input1)
            x = Dropout(0.2)(x)
        else:
            x = Dense(number_of_neurons_in_each_hlayer[i], activation= 'elu')(x)
            x = Dropout(0.2)(x)
    
    output = Dense(23, activation='sigmoid')(x)

    model = Model(inputs=[input1], outputs=output)
    return model 
    
def train_model(input1, y, x_val, y_val):

    model = get_model()

    verbose, epochs, batch_size = 1, 15, 128
    model.compile(loss='binary_crossentropy', optimizer=Adam(lr=1e-3), metrics=['accuracy'])
	# fit network
    model.fit(input1, y, epochs=epochs, batch_size=batch_size, verbose=verbose, validation_data = (x_val, y_val))
    model.save('/home/muhitemon/DeepMRG/bitscore_distribution_with_ref_exp/model_with_bitscore_distribution_against_66_exp_class_from_single_linkage40.h5')

def find_ref_exp_gene_list_and_mrg_per_class(ref_exp_fasta_file):

    ref_exp_mrg_list = []
    ref_exp_mrg_per_class = {}

    for sequence in SeqIO.parse(ref_exp_fasta_file, "fasta"):

        prot_description = str(sequence.description)
        this_exp_prot_id = prot_description.split('|')[0].strip()
        ref_exp_mrg_list.append(this_exp_prot_id)

        this_prot_label = prot_description.split('|')[-1]
        if this_prot_label.lower() == "multimetal":
            this_prot_label = prot_description.split('|')[-2]

        this_prot_class = find_ref_exp_mrg_class(this_prot_label)

        if this_prot_class not in ref_exp_mrg_per_class:
            ref_exp_mrg_per_class[this_prot_class] = [this_exp_prot_id]
        else:
            ref_exp_mrg_per_class[this_prot_class].append(this_exp_prot_id)
        
    return ref_exp_mrg_list, ref_exp_mrg_per_class # list, dict 

def compute_bitscore_feature_mat_and_label(bitscore_distribution, seq_fasta_file, ref_exp_mrg_list, ref_exp_mrg_per_class_dict):

    x = []
    y = []
    bitscore_if_no_alignment_found = 5 # defalut bitscore is set to 5 if no alignment report is found

    for sequence in SeqIO.parse(seq_fasta_file, "fasta"):

        prot_description = str(sequence.description)
        this_prot_id = prot_description.split('|')[0].strip()
        
        # something to do here to cover the case of proteins having no alignment at all from DIAMOND

        this_prot_bitscore_distribution = bitscore_distribution[this_prot_id] # this_prot_bitscore_distribution is a dict where key = ref_exp_gid and value = corresponding bitscore against the ref_exp_gid
        this_prot_x = [] # this refers to current protein's bitscore distribution against 485 ref exp MRGs

        sorted_ref_exp_mrg_classes = sorted(ref_exp_mrg_per_class_dict.keys())

        for i in sorted_ref_exp_mrg_classes:
            this_class_ref_exp_mrgs = ref_exp_mrg_per_class_dict[i]
            l = len(this_class_ref_exp_mrgs)
            sum_bitscore_of_this_class = 0
            for ref_exp_mrg_id in this_class_ref_exp_mrgs:
                if ref_exp_mrg_id in this_prot_bitscore_distribution:
                    sum_bitscore_of_this_class+=this_prot_bitscore_distribution[ref_exp_mrg_id]
                else:
                    sum_bitscore_of_this_class+=bitscore_if_no_alignment_found
            avg_bitscore_of_this_class = sum_bitscore_of_this_class / l
            this_prot_x.append(avg_bitscore_of_this_class)


        this_prot_label = prot_description.split('|')[-1]
        if this_prot_label.lower() == "multimetal":
            this_prot_label = prot_description.split('|')[-2]

        x.append(this_prot_x)
        y.append(create_lablels_of_protein_seq(this_prot_label))

    x = np.array(x)
    y = np.array(y)

    return x, y
    
if __name__ == "__main__":
    

    ref_seq_path = '/home/muhitemon/DeepMRG/processed_MRG-DB/MRG_DB_metal_wise_exp.fasta'
    ref_exp_mrg_list, ref_exp_mrg_per_class_dict = find_ref_exp_gene_list_and_mrg_per_class(ref_seq_path)
    
    
    train_seq_path = "/home/muhitemon/DeepMRG/full_len_prot_train/train_full_len_pred_per_class_single_linkage40.fasta"
    train_bitscore_distribution_file = '/home/muhitemon/DeepMRG/bitscore_distribution_with_ref_exp/train_bitscore_distribution_with_ref_exp_from_single_linkage40.json'

    with open(train_bitscore_distribution_file) as json_file:
        train_bitscore_distribution = json.load(json_file)

    #print(len(train_bitscore_distribution))

    val_seq_path = "/home/muhitemon/DeepMRG/full_len_prot_validation/validation_full_len_pred_per_class_single_linkage40.fasta"
    val_bitscore_distribution_file = '/home/muhitemon/DeepMRG/bitscore_distribution_with_ref_exp/validation_bitscore_distribution_with_ref_exp_from_single_linkage40.json'
    
    with open(val_bitscore_distribution_file) as json_file:
        val_bitscore_distribution = json.load(json_file)
    
    #print(len(val_bitscore_distribution))

    train_x, train_y = compute_bitscore_feature_mat_and_label(train_bitscore_distribution, train_seq_path, ref_exp_mrg_list, ref_exp_mrg_per_class_dict)
    val_x, val_y = compute_bitscore_feature_mat_and_label(val_bitscore_distribution, val_seq_path, ref_exp_mrg_list, ref_exp_mrg_per_class_dict)
    
    #print(train_x.shape)
    #print(val_x.shape)

    train_x_normalized = minmax_scale(train_x, axis=1) # sample (row) wise normalization instead of column wise
    val_x_normalized = minmax_scale(val_x, axis=1) # sample (row) wise normalization instead of column wise

    train_model(train_x_normalized, train_y, val_x_normalized, val_y)