# This script will run the CA updating procedure after each BRER model is created
from updating_utils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--latest-structure', type=str,
                    help='Latest BRER structure')
parser.add_argument('-r', '--run-index', type=int,
                    help='The ensemble number for each group, to separate files for each ensemble')
args = parser.parse_args()

#enter in the relevant data here
latest_structure = get_last_model(args.latest_structure)[0]  #latest structure from BRR
label_pairs = ['148_266', '66_148']    #label pairs in DEER experiments. May need to be int in functions
ca_dist_filename = f'ca_dist_dictionary_{args.run_index}'
ca_index_filename = 'ca_index_dictionary'
json_filename = f'pair_data/pair_data_{args.run_index}.json'
exp_data = glob.glob('b2*test.txt')    #get array of experimental data traces
distr_bin = f'modelled_ntx_bin_{args.run_index}' # prefix name of file in which to hold modelled distributions
learn_rate = 0.2

#process latest structure and get updated CA distance
for exp in exp_data:
    for label_pair in label_pairs:
        if label_pair in exp:
            model_ntx_update_ca(latest_structure, label='I1M', label_pair=label_pair, 
                                exp_data=exp, ca_bin=ca_dist_filename, distr_bin=distr_bin,
                                ens_num=args.run_index, learn_rate=learn_rate)
            
# create and update the pair_dist.json file
make_pair_data_file(ca_dist_filename, ca_index_filename, json_filename)