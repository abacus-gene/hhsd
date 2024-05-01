'''
EVALUATION OF AND IMPLEMENTATION OF 
PROPOSED CHANGES TO THE SPECIES DELIMITAITON
'''



import pandas as pd
import numpy as np
from typing import Dict, Literal

from .customtypehints import AlgoMode, CfileParam, gdi, NodeName, MigrationRates, Bound
from .module_ete3 import Tree, TreeNode
from .module_tree import get_node_pairs_to_modify, get_attribute_filtered_tree, get_current_leaf_species, get_iteration, add_attribute_gdi, get_attribute_filtered_imap, tree_to_newick
from .module_helper import flatten
from .module_migration import check_migration_reciprocal
from .module_gdi_numeric import get_pg1a_numerical
from .module_gdi_simulate import get_pg1a_from_sim
from .module_msa_imap import imapfile_write
from .module_parameters import ParamGDI


def get_gdi_values(
        tree:           Tree, 
        mode:           AlgoMode,
        migration_df:   MigrationRates,
        bound:          Bound
        ) ->            Dict[NodeName, gdi]:

    '''
    Get gdi values for nodes based on the parameters of the MSC(M) model

    tree is the tree datastructure holding the species delimitation
    mode is the mode of the algorithm, either 'merge' or 'split'
    migration_df is the dataframe holding the migration rates
    bound is the bound of the gdi value to be calculated, either 'lower', 'mean', or 'upper'

    '''

    # get the mode pairs for which the gdi needs to be calculated
    node_pairs_to_mod = get_node_pairs_to_modify(tree, mode)

    # create empty list of gdi values for the relevant nodes
    gdi_values = {node.name:None for node in flatten(node_pairs_to_mod)}

    # iterate through the node pairs
    for pair in node_pairs_to_mod:
        # if nodes are not involved in any migration events, or only involved in reciprocal migration events, calculate the gdi numerically
        if check_migration_reciprocal(pair[0], pair[1], migration_df) == True:
            for node in pair:
                # calculate the gdi values
                gdi_values[node.name] = get_pg1a_numerical(node, migration_df, bound)

        # otherwise, use simulation to calculate the gdi
        else:
            for node in pair:
                # calculate the gdi values
                gdi_values[node.name] = get_pg1a_from_sim(node, tree, mode, migration_df, bound)

    return gdi_values


def calculate_gdi(
        tree:           Tree, 
        mode:           AlgoMode,
        migration_df:   MigrationRates
        ) ->            Tree: 

    '''
    Calculate gdi values and add them as attributes to the Tree
    '''

    print()

    gdi_values_lower = get_gdi_values(tree, mode, migration_df, bound='lower')
    gdi_values_mean  = get_gdi_values(tree, mode, migration_df, bound='mean')
    gdi_values_upper = get_gdi_values(tree, mode, migration_df, bound='upper')

    # This print statement is used to clear the line of the previous print statements detailing the simulaltion work
    print("                                                                                                                 ", end="\r")

    # Combine the results into a single ParamGDI object
    gdi_values = {node: ParamGDI(gdi_low=gdi_values_lower[node], gdi_mean=gdi_values_mean[node], gdi_high=gdi_values_upper[node]) for node in gdi_values_lower.keys()}

    tree = add_attribute_gdi(tree, mode, gdi_values)

    return tree




# DECIDE WHETER TO ACCEPT OR REJECT PROPOSAL
def node_pair_decision(
        node_1:         TreeNode,
        node_2:         TreeNode,
        cf_dict:        CfileParam
        ) ->            None: # in-place modification of node attributes

    '''
    Decide to accept or reject a merge/split proposal
    '''

    mode:AlgoMode = cf_dict['mode']
    gdi_threshold = cf_dict['gdi_threshold']

    if   mode == 'merge':
        # if at least one of the 2 populations has low gdi indicating non-species status, merge the nodes
        if gdi_threshold == None or (node_1.gdi_mean <= gdi_threshold) or (node_2.gdi_mean <= gdi_threshold): 
            node_1.species = False; node_1.modified = True
            node_2.species = False; node_2.modified = True
        else:
            node_1.modified = False; node_2.modified = False


    elif mode == 'split':
        # if both populations have high gdi indiciating, disctinct species status, split them into 2 species
        if gdi_threshold == None or ((node_1.gdi_mean >= gdi_threshold) and (node_2.gdi_mean >= 0.5)) or ((node_2.gdi_mean >= gdi_threshold) and (node_1.gdi_mean >= 0.5)):
            node_1.species = True;  node_1.modified = True
            node_2.species = True;  node_2.modified = True
        else:
            node_1.modified = False; node_2.modified = False



def print_decision_feedback(
        node_pairs_to_modify,
        tree:                   Tree,
        cf_dict:                CfileParam
        ) ->                    None: # prints to screen, and writes files

    '''
    - Prints feedback about each of the node pairs where modifications where proposed.\\
    - Writes .csv files containing the names of node pairs and their gdi scores.       
    '''

    # function to collect parameters about the node pair relevant to the decision
    def node_pair_feedback(node_1, node_2, cf_dict):
        results = {
            "node 1":node_1.name,
            "GDI 1": node_1.gdi_mean,
            "lower bound 1" : node_1.gdi_low,
            "upper bound 1" : node_1.gdi_high,
            "node 2":node_2.name,
            "GDI 2": node_2.gdi_mean,
            "lower bound 2" : node_2.gdi_low,
            "upper bound 2" : node_2.gdi_high,
            
            f"{cf_dict['mode']} accepted?": node_1.modified}

        return results

    # collect feedback
    feedback = []
    for pair in node_pairs_to_modify:
        feedback.append(node_pair_feedback(pair[0], pair[1], cf_dict))

    df = pd.DataFrame.from_dict(feedback)
    df.rename(columns={'lower bound 1':"2.5% HPD", 'upper bound 1':"97.5% HPD",'lower bound 2':"2.5% HPD", 'upper bound 2':"97.5% HPD"}, inplace=True)
    df.to_csv("decision.csv", index=False)

    # print each node pair, the gdi, and whether the proposal as accepted
    print(f"\n> Proposal results:\n")
    print(df.to_string(index=False, max_colwidth=36, justify="start"))
    
    # print the names of currently accepted species
    print(f"\nNumber of species after iteration {get_iteration(tree)}:  {len(get_current_leaf_species(tree))}")
    print(str(get_current_leaf_species(tree))[1:-1])

    # print current topology as ASCII
    print(get_attribute_filtered_tree(tree, "species", newick=False))



## FINAL WRAPPER FUNCTION
def tree_modify_delimitation(
        tree:       Tree,
        cf_dict:    CfileParam,
        ) ->        Tree:     

    '''
    Modify the tree datastructure holding the species delimitation accoring to the inferred gdi parameters, 
    and the gdi thresholds set by the user.
    '''

    # get the node pairs where modifications are being proposed
    node_pairs_to_modify = get_node_pairs_to_modify(tree, cf_dict['mode'])

    # modify node pairs according to criteria
    for pair in node_pairs_to_modify:
        node_pair_decision(pair[0], pair[1], cf_dict)
    
    print_decision_feedback(node_pairs_to_modify, tree, cf_dict)

    # output the current imap of accepted species
    result_imap = get_attribute_filtered_imap(tree, attribute='species')
    imapfile_write(result_imap, 'RESULT_IMAP.txt')

    # output the currently accepted newick tree
    result_tree = get_attribute_filtered_tree(tree, attribute='species')
    with open('RESULT_TREE.txt', 'w') as f:
        f.write(result_tree)

    return tree