import chilife as xl
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import os
from natsort import os_sorted
import glob
import pickle

def gro2pdb(structure):
    '''
    This function creates a pdb from a gro file

    This may not be needed at all
    '''

def initialize_files(starting_model, label_pairs, ca_dist_filename='ca_dist_dict', ca_index_filename='ca_index_dict'):
    '''
    Initialize the files that will hold the C-alphas distances and indices. The keys for these two dictionaries should
    be identical for use later on.
    '''
    ca_dist_dict = {}
    ca_index_dict = {}
    for pair in label_pairs:
        site1, site2 = pair.split('_')
        u=mda.Universe(starting_model)
        ca_1 = u.select_atoms(f'resid {site1} and name CA')
        ca_2 = u.select_atoms(f'resid {site2} and name CA')
        res1, res2, ca = dist(ca_1, ca_2)

        ca_dist_dict[f'{pair}'] = []
        ca_dist_dict[f'{pair}'].append(ca)

        ca_index_dict[f'{pair}'] = []
        ca_index_dict[f'{pair}'].append([ca_1.indices[0]+1, ca_2.indices[0]+1])
    
    with open(f'{ca_dist_filename}.pickle', 'wb') as handle:
        pickle.dump(ca_dist_dict, handle)
        handle.close()

    with open(f'{ca_index_filename}.pickle', 'wb') as handle:
        pickle.dump(ca_index_dict, handle)
        handle.close()

def make_pair_data_file(ca_dictionary, ca_indices, json_file_name='pair_data.json'):
    '''
    Given a C-alpha dictionary and indices from the starting structure, this utility will create the appropriately
    formatted pair_data.json file for input to BRER. 
    '''
    with open(f'{ca_dictionary}', 'rb') as file:    #unpickling C-alpha distance dictionary
        ca_dictionary = pickle.load(file)
    with open(f'{ca_indices}', 'rb') as file:       #unpickling C-alpha index dictionary
        ca_indices = pickle.load(file)
    
    keys = list(ca_dictionary.keys())

    json_template = ['{\n']
    json_file = open("{}".format(json_file_name), "w")
    json_file.writelines(json_template)
    json_file.close()
    count = 0

    for key in keys:
        ca1_index, ca2_index = ca_indices[key]
        ca_dist = ca_dictionary[key][-1] #pulls the latest CA distance update
        if count < len(keys):
            name = f'{key}'
            sites = f'{ca1_index}, {ca2_index}'
            dist = ca_dist/10 #this needs to be in nm for BRER/Gromacs
            dist = dist.tolist()
            prob = 1            # 100% probability to pick this point (it's the only point lol)
            prob=prob.tolist()
            
            json_template = [
            '   "{}":'.format(name),' {\n',
            '       "sites": [{}],\n'.format(sites), 
            '       "name": "{}",\n'.format(name),
            '       "distribution": {},\n'.format(prob),
            '       "bins": {}\n'.format(dist),
            '    },\n']
            
            json_file = open("{}".format(json_file_name), "a")
            json_file.writelines(json_template)
            json_file.close()
            
        elif count == len(keys):
            name = f'{key}'
            sites = f'{ca1_index}, {ca2_index}'
            dist = ca_dist/10 #this needs to be in nm for BRER/Gromacs
            dist = dist.tolist()
            prob = 1            # 100% probability to pick this point (it's the only point lol)
            prob=prob.tolist()
            
            json_template = [
            '   "{}":'.format(name),' {\n',
            '       "sites": [{}],\n'.format(sites), 
            '       "name": "{}",\n'.format(name),
            '       "distribution": {},\n'.format(prob),
            '       "bins": {}\n'.format(dist),
            '    }\n',
            '}']
            
            json_file = open("{}".format(json_file_name), "a")
            json_file.writelines(json_template)
            json_file.close()




def get_last_model(directory):
    '''
    Searches and locates the last model produced by BRER. This will couple with argparse so that the current BRER 
    folder will be loaded and the members of that folder will then be searched, gathered, and ordered with natsort.
    Then, the last file in the string will be chosen, which is the last model created, and labeled.

    dir: string, directory to search
    '''
    _dir = os_sorted(glob.glob(directory+'/[!state.json]*')) #element the json file from the search
    _dir = os.path.join(_dir[-1], 'production/*.gro')         #take last file in dir and join with string to get last model
    return _dir 

def model_ntx_update_ca(structure, label, label_pair, exp_data, distr_bin, ca_bin='ca_bin', learn_rate=0.1, **kwargs):
    '''
    This function creates a modelled nitroxide distribution from and input structure,
    specifically from BRER. It then Updates the pair_data file with new C-alpha distances in which to bias the protein. 
    This is based on the positive residuals from exp - ntxd_model_data. The weighted average distance of the poistive 
    residuals is computed, as well as the weighted average of the cumulative model_data. The difference of these two values
    are multipled by a learning rate and posed as the next bias distance.
    
    structure: string, path to structure/model in PDB/GRO format
    
    label: string, i.e. 'I1M' or 'R1M'
    site1/2 = int, residue to label

    distr_bin = str, file name in which the modelled nitroxide distributions are held
    
    r = array, specifies x-axis of label
    '''
    u=mda.Universe(structure)
    exp_data = np.loadtxt(exp_data) #data needs to be loaded as array with first vector as r values
    r = exp_data[0] #r values need to be identical between experiment and modelled nitroxide for this to be accurate
    distr_bin = distr_bin+'_'+label_pair+'_'+label+'.txt'   #name for distribution bin of a given label pair
    with open(f'{ca_bin}', 'rb') as file:    #unpickling C-alpha distance dictionary
        ca_dictionary = pickle.load(file) #load the ca distance dictionary

    site1, site2 = label_pair.split('_')
    
    SL1 = xl.SpinLabel(protein=u, label=label, site=site1)
    SL2 = xl.SpinLabel(protein=u, label=label, site=site2)

    #need to add in error handling for the spin label, especially V1X
    #also need to create a function to do this for each label pair. Might not have to be part of utils
    traj, de = xl.repack(u, SL1, SL2,
                            repetitions=2500,
                            temp=295,
                            off_rotamer=False,
                            repack_radius=10) 
    SL1r = xl.SpinLabel.from_trajectory(traj, site=int(site1), burn_in=1000, spin_atoms=SL1.spin_atoms)
    SL2r = xl.SpinLabel.from_trajectory(traj, site=int(site2), burn_in=1000, spin_atoms=SL2.spin_atoms)

    P = np.array(xl.distance_distribution(SL1r, SL2r, r))
    
    if os.path.exists(distr_bin) == False:
         np.savetxt(distr_bin, P)
         updated = P
    else:
         d = np.loadtxt(distr_bin)
         updated = np.vstack(d, P)
         np.savetxt(distr_bin, updated)

    #Update C-alpha distance 
    if len(updated.shape) == 1:
        residual = exp_data-updated
        residual[residual<0]=0  #select for positive residuals

        #if updated is only one vector long, then the initial CA distance has not been measured yet
        ca_1 = u.select_atoms(f'resid {site1} and name CA')
        ca_2 = u.select_atoms(f'resid {site2} and name CA')
        res1, res2, ca = dist(ca_1, ca_2)

    else:
        updated = updated.sum(axis=0)
        residual = exp_data-updated
        residual[residual<0]=0  #select for positive residuals

    res_w_avg = np.average(r, weights=residual/sum(residual))
    mod_w_avg = np.average(r, weights=updated/sum(updated))

    prev_ca = ca_dictionary[label_pair][-1]
    new_ca = prev_ca + (learn_rate*(res_w_avg-mod_w_avg)) 

    #append new ca distance into CA dictionary
    ca_dictionary[label_pair].append(new_ca)
    with open(f'{ca_bin}', 'wb') as file:
        pickle.dump(ca_dictionary, file)
        file.close()




def init_dist(starting_model, exp_data):
    '''
    Get the initial C-alpha distance to run for BRER. This will measure the starting structure CA distance, model nitroxide,
    and provide a new CA distance as the first run from the model_ntx_update_ca function. 
    '''
    u=mda.Universe(starting_model)
    exp_data = np.loadtxt(exp_data) #data needs to be loaded as array with first vector as r values
    r = exp_data[0] #r values need to be identical between experiment and modelled nitroxide for this to be accurate
    
    SL1 = xl.SpinLabel(protein=u, label=label, site=site1)
    SL2 = xl.SpinLabel(protein=u, label=label, site=site2)

    #need to add in error handling for the spin label, especially V1X
    #also need to create a function to do this for each label pair. Might not have to be part of utils
    traj, de = xl.repack(u, SL1, SL2,
                            repetitions=2500,
                            temp=295,
                            off_rotamer=False,
                            repack_radius=10) 
    SL1r = xl.SpinLabel.from_trajectory(traj, site=int(site1), burn_in=1000, spin_atoms=SL1.spin_atoms)
    SL2r = xl.SpinLabel.from_trajectory(traj, site=int(site2), burn_in=1000, spin_atoms=SL2.spin_atoms)

    P = np.array(xl.distance_distribution(SL1r, SL2r, r))
    
    if os.path.exists(distr_bin) == False:
         np.savetxt(distr_bin, P)
    else:
         d = np.loadtxt(distr_bin)
         updated = np.vstack(d, P)
         np.savetxt(distr_bin, updated)
         
def update_ca(exp_data, ntxd_model_data, pair_data):
    '''
    Updates the pair_data file with new C-alpha distances in which to bias the protein. This is based on the positive
    residuals from exp - ntxd_model_data. The weighted average distance of the poistive residuals is computed, as well 
    as the weighted average of the cumulative model_data. The difference of these two values are multipled by a learning
    rate and added to a randomly drawn C


    '''
    exp_data = np.loadtxt(exp_data)