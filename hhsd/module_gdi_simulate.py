'''
INFER GDI VALUES BY SIMULATING GENETREES UNDER THE MSC+M MODEL

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

import re
import copy
import subprocess
import os

from .customtypehints import BppCfile, BppCfileParam, GeneTrees, AlgoMode, MigrationRates, NodeName
from .module_ete3 import Tree, TreeNode
from .module_helper import readlines, dict_merge, get_bundled_bpp_path
from .module_bpp import bppcfile_write
from .module_tree import get_attribute_filtered_tree, add_attribute_tau_theta, ensure_taus_valid
from .module_bpp_readres import MSCNumericParamEstimates, NumericParam


def tree_to_extended_newick(
        tree:   Tree,
        ) ->    str:

    '''
    'tree' is an ete3 Tree object which contains the topology, and the tau and theta values as node attributes.
    The output is a newick tree that also contains information about the tau and theta values at each node.
    This newick tree is used as input to bpp --simulate
    '''

    root = tree.get_tree_root()

    # create the extended newick version of the tree topology which contains the tau and theta values
    tree_str = tree.write(features = ['tau', 'theta'],format=1)

    # filter out extra material not related to required parameters
    tree_str = re.sub(r':1\[&&NHX', '', tree_str)
    tree_str = re.sub(f':theta=', ' #', tree_str)
    tree_str = re.sub(f':tau=None', '', tree_str)
    tree_str = re.sub(f':tau=', ' :', tree_str)
    tree_str = re.sub(r'\]', '', tree_str)
    tree_str = re.sub(r'\)', ') ', tree_str)
    tree_str = re.sub(r'\(', ' (', tree_str)
    
    # add in data corresponding to root node, which is not added in by ete3 for some reason
    tree_str = re.sub(';', f'{root.name} :{root.tau} #{root.theta};', tree_str)

    return tree_str

def get_migration_events(
        migration_df: MigrationRates,                  
        ) -> str:

    """
    Get the migration events and rates from the migration dataframe to append to the control file.

    The migration events are appended to the hhsd control file in the following format:
    'migration = n
    source destination W
    """    

    mig_events = f'migration = {len(migration_df["source"])}\n {migration_df.to_string(header = False, index = False)}'

    return mig_events


# default parameters for a 'bpp --simulate' control file used for simulating gene trees
default_BPP_simctl_dict:BppCfileParam = {
    'seed':                 '1111',
    'treefile':             'MyTree.tre', 
    'Imapfile':             'MyImap.txt', 
    'species&tree':         None, 
    'popsizes':             None, 
    'newick':               None,
    'loci&length':          '1000 50',
}


# create the bpp --simulate cfile for simulating gene trees
def create_simulate_cfile(
        node:           TreeNode,
        tree:           Tree, 
        mode:           AlgoMode, 
        migration_df:   MigrationRates,
        ) ->            None: # writes control file to disk

    '''
    - 'tree' is an ete3 Tree object.
    - 'mode' specifies whether the algo is running in merge or split mode.
    - 'migration_df' is the DataFrame object containing the source, destination, and rate (M) for all migration events.
    - 'bound' is the bound of the gdi value to be calculated, either 'lower', 'mean', or 'upper'.

    the function writes a 'bpp --simulate' control file to disk specifying the parmeters of the simulation. 
    All populations in the simulation generate two sequences, as this facilitates the estiamtion of the gdi from gene trees 
    (performed in 'get_gdi_from_sim').
    '''

    # get tree object needed to create simulation 
    sim_tree = get_attribute_filtered_tree(tree, mode, newick=False)
    # ensure descendants are younger than ancestors, and modify tree if needed
    sim_tree = ensure_taus_valid(sim_tree)
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
    mig_events = get_migration_events(migration_df)
    with open("sim_ctl.ctl", "a") as myfile: 
        myfile.write(mig_events)



def run_BPP_simulate(
        control_file:   BppCfile,  
        ) ->            None: # handles the bpp subprocess
    
    '''
    Use 'bpp --simulate' to sample gene trees from a given MSC+M model
    '''

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
            break


# final wrapper function to simulation gene trees according to the given MSC+M model
def genetree_simulation(
        node:           TreeNode,
        tree:           Tree, 
        mode:           AlgoMode, 
        migration_df:   MigrationRates,
        ) ->            GeneTrees: 

    '''
    Handle the file system operations, and bpp control file creation to simulate gene trees. Return the gene trees as a list
    '''

    # create temporary directory to store bpp --simulate output
    os.mkdir('genetree_simulate')
    os.chdir('genetree_simulate')

    # write the cfile to disk
    create_simulate_cfile(node, tree, mode, migration_df)
    
    # run bpp --simulate
    run_BPP_simulate('sim_ctl.ctl')
    
    # read the gene trees from the output file
    try:
        all_genetrees = readlines('MyTree.tre')
    except:
        raise ValueError("Error in simulating gene trees. Please check the /genetree_simulate folder for more information.")

    os.remove('MyTree.tre')
    os.remove('MyImap.txt')
    os.remove('sim_ctl.ctl')
    os.chdir('..')
    os.rmdir('genetree_simulate')

    return all_genetrees


def pg1a_from_genetrees(
        node:           TreeNode,
        tau_AB:         float, 
        all_genetrees:  GeneTrees
        ) ->            float:
    
    '''
    Get P(G1A) of a given TreeNode from the simulated data.

    P(G1A) is probability of the topology ((a1, a2), b1) before the populations split.

    This definition allows us to estimate P(G1) from simulated gene tree topologies. After simulating many genetrees for a 
    fully specified MSC+M model with two sequences from population A, and one from the sister population B, P(G1A) of A can be estimated by 
    counting the proportion of gene trees where the topology ((a1, a2), b1) is observed before the populations split.
    '''

    node_name = str(node.name)
    seq_name = node_name.lower()


    # create the regex corresponding to the required topology
    correct_genetree = f'\({seq_name}[12]\^{node_name}[:][0][.][\d]#,{seq_name}[12]\^{node_name}:[0][.][\d]#\)'
    correct_genetree = re.sub('#', '{6}', correct_genetree) # this is needed due to no '{''}' characters being allowed within f strings
    
    # find all occurrences of the correct topology
    g1_genetrees = [re.search(correct_genetree, element) for element in all_genetrees]
    g1_genetrees = [element.group(0) for element in g1_genetrees if element]

    # isolate the time at which each correct topology is achieved
    num_match = (re.search('0.\d{6}', g1_genetrees[0])); num_start = num_match.span()[0]; num_end = num_match.span()[1]
    times = [float(element[num_start:num_end]) for element in g1_genetrees]

    # isolate the ocurrences that are after the split time for the populations
    before_split = [element for element in times if element < tau_AB]

    # gdi is the proportion of the loci where this topology is observed before the populations split
    pg1a = len(before_split)/len(all_genetrees)

    return pg1a

def get_pg1a_from_sim(
        node:           TreeNode,
        tree:           Tree,
        mode:           AlgoMode,
        numeric_param:  MSCNumericParamEstimates,
        ) ->            NumericParam:
    
    '''
    Get P(G1A) of a given TreeNode by simulating trees and counting the proportion of trees with the correct topology.
    Perform the 1000 replicate simulations needed to establish a sample from the distribution over the gdi. Various
    statistics [mean, confidence intervals, etc] can then be estimated from this sample.
    '''

    ancestor_node:NodeName = str(node.up.name)

    results = []
    # run the 1000 replicate gdi estimations
    for i in range(1000):
        print(f"inferring gdi for '{node.name}' using gene tree simulation ({i+1}/1000)...                        ", end="\r")

        # sample the mcmc values
        tau_dict        = numeric_param.sample_tau(i)
        theta_dict      = numeric_param.sample_theta(i)
        migration_df    = numeric_param.sample_migparam(i)

        # create a new tree object with tau and theta values corresponding to the newly sampled values
        tree_copy = copy.deepcopy(tree)
        tree_copy = add_attribute_tau_theta(tree_copy, tau_dict, theta_dict) 

        # simulate the gene trees
        genetrees = genetree_simulation(node, tree_copy, mode, migration_df)

        # get the time at which the populations split
        tau_AB = tau_dict[ancestor_node]

        # get P(G1A)
        results.append(pg1a_from_genetrees(node, tau_AB, genetrees))

    print("                                                                                                       ", end='\r')

    return NumericParam(results)