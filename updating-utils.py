import chilife as xl
import numpy as np
import MDAnalysis as mda
import os
from natsort import os_sorted
import glob

def gro2pdb(structure):
    '''
    This function creates a pdb from a gro file

    This may not be needed at all
    '''


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


def model_ntx_update_ca(structure, label, site1, site2, exp_data, distr_bin, **kwargs):
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
    
    SL1 = xl.SpinLabel(protein=u, label=label, site=site1)
    SL2 = xl.SpinLabel(protein=u, label=label, site=site2)

    #need to add in error handling for the spin label, especially V1X
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

    #Update C-alpha distance 

def init_dist(starting_model, exp_data):
    '''
    Get the initial C-alpha distance to run for BRER. This will measure the starting structure CA distance, model nitroxide,
    and provide a new CA distance as the first run from the model_ntx_update_ca function. 
    '''
         
def update_ca(exp_data, ntxd_model_data, pair_data):
    '''
    Updates the pair_data file with new C-alpha distances in which to bias the protein. This is based on the positive
    residuals from exp - ntxd_model_data. The weighted average distance of the poistive residuals is computed, as well 
    as the weighted average of the cumulative model_data. The difference of these two values are multipled by a learning
    rate and added to a randomly drawn C


    '''
    exp_data = np.loadtxt(exp_data)