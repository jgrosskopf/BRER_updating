# This script will run the CA updating procedure after each BRER model is created
from updating_utils import *
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--latest-structure', type=str,
                    help='Latest BRER structure')
args = parser.parse_args()

#enter in the relevant data here
latest_structure = get_last_model(args.latest_structure)  #latest structure from BRR
label_pairs = ['148_266', '148_228']    #label pairs in DEER experiments. May need to be int in functions
ca_dist_filename = 'ca_dist_dictionary'
ca_index_filename = 'ca_index_dictionary'
json_filename = 'pair_data.json'
exp_data = glob.glob('expdata*.txt')    #get array of experimental data traces
distr_bin = 'modelled_ntx_bin' # prefix name of file in which to hold modelled distributions
learn_rate = 0.1

#process latest structure and get updated CA distance
for exp in exp_data:
    for label_pair in label_pairs:
        if label_pair in exp:
            new_ca = model_ntx_update_ca(latest_structure, label='I1M', label_pair=label_pair, 
                                         exp_data=exp, ca_bin=ca_dist_filename, learn_rate=learn_rate)
            
# create and update the pair_dist.json file
make_pair_data_file(ca_dist_filename, ca_index_filename, json_filename)