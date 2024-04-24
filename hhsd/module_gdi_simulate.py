'''
INFER GDI VALUES BY SIMULATING GENETREES UNDER THE MSC+M MODEL
'''

import re
import copy
import subprocess
import os

import numpy as np

from .classes import BppCfile, BppCfileParam, GeneTrees, gdi, AlgoMode, MigrationRates
from .module_ete3 import Tree, TreeNode
from .module_helper import readlines, dict_merge, get_bundled_bpp_path
from .module_bpp import bppcfile_write
from .module_tree import get_attribute_filtered_tree


## INFERENCE OF GDI FROM GENETREES
'''
The functions in this section are responsible for outputing the list of 10^6 gene tree topologies with 
associated branch lengths. 

A fully specified MSC+M model consists of the following:

    1) Tree topology
    2) Branch lengths (tau)
    3) Effective population sizes (theta)
    4) Source and destination of migration events
    5) Rate corresponding to each migration event (M)

Such a model defines the joint distribution of gene tree topolgies and coalescence times,
and simulation can be used to sample this distribution. 
'''

def tree_to_extended_newick(
        tree:   Tree
        ) ->    str:

    '''
    'tree' is an ete3 Tree object which contains the topology, and the tau and theta values as node attributes.
    The output is an extended newick tree that also contains information about the tau and theta values.
    '''

    # create the extended newick version of the tree topology which contains the tau and theta values
    tree_str = tree.write(features = ['tau', 'theta'],format=1)

        # filter out extra material not related to required parameters
    tree_str = re.sub(r':1\[&&NHX', '', tree_str)
    tree_str = re.sub(':theta=', ' #', tree_str)
    tree_str = re.sub(':tau=None', '', tree_str)
    tree_str = re.sub(':tau=', ' :', tree_str)
    tree_str = re.sub(r'\]', '', tree_str)
    tree_str = re.sub(r'\)', ') ', tree_str)
    tree_str = re.sub(r'\(', ' (', tree_str)
        # add in data corresponding to root node, which is not added in by ete3 for some reason
    root = tree.get_tree_root()
    tree_str = re.sub(';', f'{root.name} :{root.tau} #{root.theta};', tree_str)

    return tree_str


# default parameters for a 'bpp --simulate' control file used for simulating gene trees
default_BPP_simctl_dict:BppCfileParam = {
    'seed':                 '1111',
    'treefile':             'MyTree.tre', 
    'Imapfile':             'MyImap.txt', 
    'species&tree':         None, 
    'popsizes':             None, 
    'newick':               None,
    'loci&length':          '1000000 500',
}


# create the bpp --simulate cfile for simulating gene trees
def create_simulate_cfile(
        node:           TreeNode,
        tree:           Tree, 
        mode:           AlgoMode, 
        migration_df:   MigrationRates
        ) ->            None: # writes control file to disk

    '''
    - 'tree' is an ete3 Tree object.
    - 'mode' specifies whether the algo is running in merge or split mode.
    - 'migration_df' is the DataFrame object containing the source, destination, and rate (M) for all migration events.

    the function writes a 'bpp --simulate' control file to disk specifying the parmeters of the simulation. 
    All populations in the simulation generate two sequences, as this facilitates the estiamtion of the gdi from gene trees 
    (performed in 'get_gdi_from_sim').
    '''

    # get tree object needed to create simulation 
    sim_tree = get_attribute_filtered_tree(tree, mode, newick=False)
    leaf_names = set([leaf.name for leaf in sim_tree])
    node_name = node.name
    sister_name = node.get_sisters()[0].name

    # infer the parameters of the simulation dict from the tree object
    sim_dict = {}
    sim_dict['species&tree'] = f'{len(leaf_names)} {" ".join(leaf_names)}'
    sim_dict['popsizes'] = '     '
    
    for nodename in leaf_names:
        # simulate two sequences from the node of interest, and one from the sister population
        if nodename == node_name:
            sim_dict['popsizes'] += '2 '
        elif nodename == sister_name:
            sim_dict['popsizes'] += '1 '
        else:
            sim_dict['popsizes'] += '0 '

    sim_dict['newick'] = tree_to_extended_newick(sim_tree)

    # write the control dict
    ctl_dict = dict_merge(copy.deepcopy(default_BPP_simctl_dict), sim_dict)
    bppcfile_write(ctl_dict,"sim_ctl.ctl")

    # append lines corresponding to migration events and rates (simulation is only required if migraiton is present in the model)
    mig_events = f'migration = {len(migration_df["source"])}\n {migration_df.to_string(header = False, index = False)}'
    with open("sim_ctl.ctl", "a") as myfile: 
        myfile.write(mig_events)



def run_BPP_simulate(
        control_file:   BppCfile,  
        ) ->            None: # handles the bpp subprocess
    
    '''
    Use 'bpp --simulate' to sample gene trees from a given MSC+M model
    '''

    print(f"\ninferring gdi using gene tree simulation...", end="\r")

    # runs BPP in a dedicated subprocess
    process = subprocess.Popen(
        f"{get_bundled_bpp_path()} --simulate {control_file}", 
        shell = True, 
        bufsize = 1,
        stdout = subprocess.PIPE, 
        stderr = subprocess.STDOUT,
        encoding = 'utf-8', 
        errors = 'replace' 
        )

    # this is necessary so that the program does not hang while the simulations are completing
    while True:
        realtime_output = process.stdout.readline()

        # exit if process is stopped
        if realtime_output == '' and process.poll() is not None:
            print("                                                                     ", end="\r")
            break


# final wrapper function to simulation gene trees according to the given MSC+M model
def genetree_simulation(
        node:           TreeNode,
        tree:           Tree, 
        mode:           AlgoMode, 
        migration_df:   MigrationRates
        ) ->            GeneTrees: 

    '''
    Handle the file system operations, and bpp control file creation to simulate gene trees. Return the gene trees as a list
    '''

    
    os.mkdir('genetree_simulate')
    os.chdir('genetree_simulate')

    # write the cfile to disk
    create_simulate_cfile(node, tree, mode, migration_df)
    
    # run bpp --simulate
    run_BPP_simulate('sim_ctl.ctl')
    
    # read the gene trees from the output file
    all_genetrees = readlines('MyTree.tre')

    os.remove('MyTree.tre') # delete the tree file (it is very large)
    os.chdir('..')
    os.rmdir('genetree_simulate')

    return all_genetrees


def get_pg1_from_sim(
        node:           TreeNode, 
        all_genetrees:  GeneTrees
        ) ->            float:
    
    '''
    Get P(G1) of a given TreeNode from the simulated data.

    P(G1) is probability of the topology ((a1, a2), b1);

    This definition allows us to estimate P(G1) from simulated gene tree topologies. After simulating many genetrees for a 
    fully specified MSC+M model with two sequences from population A, and one from the sister population B, P(G1) of A can be estimated by 
    counting the proportion of gene trees where the topology ((a1, a2), b1) is observed.
    '''

    node_name = str(node.name)
    seq_name = node_name.lower()

    # create the regex corresponding to the required topology
    correct_genetree = f'\({seq_name}[12]\^{node_name}[:][0][.][\d]#,{seq_name}[12]\^{node_name}:[0][.][\d]#\)'
    correct_genetree = re.sub('#', '{6}', correct_genetree) # this is needed due to no '{''}' characters being allowed within f strings
    
    # find all occurrences of the correct topology
    g1_genetrees = [re.search(correct_genetree, element) for element in all_genetrees]
    g1_genetrees = [element.group(0) for element in g1_genetrees if element]

    # gdi is the proportion of the loci where this topology is observed
    gdi = len(g1_genetrees)/len(all_genetrees)

    return np.round(gdi, 2)