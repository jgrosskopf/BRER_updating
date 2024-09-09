# This is the script that will be called before BRER runs to initialize all files
from updating_utils import *
import glob

#enter in the relevant data here
starting_structure = 'start_model.gro'  #starting structure into BRER
label_pairs = ['148_266', '148_228']    #label pairs in DEER experiments. May need to be int in functions
ca_dist_filename = 'ca_dist_dictionary'
ca_index_filename = 'ca_index_dictionary'
json_filename = 'pair_data.json'
exp_data = glob.glob('expdata*.txt')    #get array of experimental data traces
distr_bin = 'modelled_ntx_bin' # prefix name of file in which to hold modelled distributions
learn_rate = 0.1

#initializing dictionaries
initialize_files(starting_structure, label_pairs, ca_dist_filename, ca_index_filename)

#process starting structure and get updated CA distance for first run
for exp in exp_data:
    for label_pair in label_pairs:
        if label_pair in exp:
            new_ca = model_ntx_update_ca(starting_structure, label='I1M', label_pair=label_pair, 
                                         exp_data=exp, ca_bin=ca_dist_filename, learn_rate=learn_rate)

#create first pair_data.json instance
make_pair_data_file(ca_dist_filename, ca_index_filename, json_filename)