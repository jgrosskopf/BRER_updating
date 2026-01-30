# This is the script that will be called before BRER runs to initialize all files and calculate initial CA distance targets

from updating_utils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--starting-structure', type=str,
                    help='Starting BRER structure')
parser.add_argument('-r', '--run-index', type=int,
                    help='The ensemble number for each group, to separate pair_data files')
args = parser.parse_args()

#enter in the relevant data here
starting_structure = args.starting_structure  #starting structure into BRER
label_pairs = ['55_139', '55_236', '55_304', '55_311', '139_220', '139_236', '139_304', '139_311', '220_311', '236_311']
ca_dist_filename = f'ca_dist_dictionary_{args.run_index}'
ca_index_filename = 'ca_index_dictionary'
json_filename = f'pair_data/pair_data_{args.run_index}.json'
exp_data = glob.glob('*final.txt')    #get array of experimental data traces
distr_bin = f'modelled_ntx_bin_{args.run_index}' # prefix name of file in which to hold modelled distributions
learn_rate = 0.2
momentum = 1

#initializing dictionaries
initialize_files(starting_structure, label_pairs, ca_dist_filename, ca_index_filename)

#process starting structure and get updated CA distance for first run
for exp in exp_data:
    for label_pair in label_pairs:
        if label_pair in exp:
            model_ntx_update_ca(starting_structure, label='V1X_sampled', label_pair=label_pair, 
                                exp_data=exp, ca_bin=ca_dist_filename, distr_bin=distr_bin,
                                ens_num=args.run_index, learn_rate=learn_rate, momentum=momentum)

#create first pair_data.json instance
make_pair_data_file(ca_dist_filename, ca_index_filename, json_filename)