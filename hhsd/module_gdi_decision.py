'''
EVALUATION OF AND IMPLEMENTATION OF 
PROPOSED CHANGES TO THE SPECIES DELIMITAITON
'''

import pandas as pd
import numpy as np
from typing import Dict, Literal

from .customtypehints import AlgoMode, CfileParam, NodeName, MigrationRates
from .module_ete3 import Tree, TreeNode
from .module_tree import get_node_pairs_to_modify, get_attribute_filtered_tree, get_current_leaf_species, get_iteration, get_attribute_filtered_imap
from .module_helper import flatten, check_numeric
from .module_migration import check_migration_reciprocal
from .module_gdi_numeric import get_pg1a_numerical
from .module_gdi_simulate import get_pg1a_from_sim
from .module_msa_imap import imapfile_write
from .module_bpp_readres import MSCNumericParamEstimates, NumericParam


def get_gdi_values(
        tree:           Tree, 
        numeric_param:  MSCNumericParamEstimates,
        mode:           AlgoMode,
        ) ->            Dict[NodeName, NumericParam]:

    '''
    Get gdi values for nodes based on the parameters of the MSC(M) model

    tree is the tree datastructure holding the species delimitation
    numeric_param holds the results of the MCMC on the MSC model
    mode is the mode of the algorithm, either 'merge' or 'split'
    '''

    # get the mode pairs for which the gdi needs to be calculated
    node_pairs_to_mod = get_node_pairs_to_modify(tree, mode)

    # create empty list of gdi values for the relevant nodes
    gdi_values = {node.name:None for node in flatten(node_pairs_to_mod)}

    # create small migration df, only used to check reciprocity
    migdf_for_reciproc_check = numeric_param.sample_migparam(0)

    # iterate through the node pairs
    for pair in node_pairs_to_mod:
        # if nodes are not involved in any migration events, or only involved in reciprocal migration events, calculate the gdi numerically
        if check_migration_reciprocal(pair[0], pair[1], mig_pattern=migdf_for_reciproc_check) == True:
            for node in pair:
                gdi_values[node.name] = get_pg1a_numerical(node, numeric_param)

        # otherwise, use simulation to calculate the gdi
        else:
            for node in pair:
                gdi_values[node.name] = get_pg1a_from_sim(node, tree, mode, numeric_param)

    return gdi_values


# DECIDE WHETER TO ACCEPT OR REJECT PROPOSAL
def node_pair_decision(
        node_1:         TreeNode,
        node_2:         TreeNode,
        gdi_values:     Dict[NodeName, NumericParam],
        cf_dict:        CfileParam
        ) ->            None: # in-place modification of node attributes

    '''
    Decide to accept or reject a merge/split proposal
    '''

    # get the thresholds (these come in the form similar to >0.4 or <=1.0)
    gdi_thresholds = cf_dict['gdi_threshold']
    threshold_1 = f'x{gdi_thresholds[0]}' 
    threshold_2 = f'x{gdi_thresholds[1]}'

    # get the mean gdis (formatted to strings to comply with the check_numeric function's expected input type)
    mean_gdi_1 = str(np.round(gdi_values[str(node_1.name)].mean(),2))
    mean_gdi_2 = str(np.round(gdi_values[str(node_2.name)].mean(),2))

    # check if either of the two possible combinations of interpreting the gdi values and thresholds evaluate to true. 
    #       for example, in merge mode, if mean_gdi_1 = 0.9, mean_gdi_2 = 0.8, then thresholds of >0.7,>0.7 will eval to true.     
    within_thresholds = (check_numeric(mean_gdi_1, threshold_1) and check_numeric(mean_gdi_2, threshold_2)) or (check_numeric(mean_gdi_1, threshold_2) and check_numeric(mean_gdi_2, threshold_1))
    
    # modify node attribues to reflect if the values where within the thresholds.
    if   cf_dict["mode"] == "merge": # in merge mode, candidates are stripped of their species status
        if within_thresholds: 
            node_1.species = False; node_1.modified = True
            node_2.species = False; node_2.modified = True
        else:
            node_1.modified = False; node_2.modified = False
    
    elif cf_dict["mode"] == "split": # in split mode, candidates are granted species status
        if within_thresholds:
            node_1.species = True;  node_1.modified = True
            node_2.species = True;  node_2.modified = True
        else:
            node_1.modified = False; node_2.modified = False
 

def print_decision_feedback(
        node_pairs_to_modify,
        tree:                   Tree,
        gdi_values:             Dict[NodeName, NumericParam],
        cf_dict:                CfileParam
        ) ->                    None: # prints to screen, and writes files

    '''
    - Prints feedback about each of the node pairs where modifications where proposed.\\
    - Writes .csv files containing the names of node pairs and their gdi scores.       
    '''

    # function to collect parameters about the node pair relevant to the decision
    def node_pair_feedback(node_1, node_2, gdi_values, cf_dict):
        node_1_name = str(node_1.name)
        node_1_gdi_mean = gdi_values[node_1_name].mean()
        node_1_gdi_hpd = gdi_values[node_1_name].hpd_bound(0.95)
        
        node_2_name = str(node_2.name)
        node_2_gdi_mean = gdi_values[node_2_name].mean()
        node_2_gdi_hpd = gdi_values[node_2_name].hpd_bound(0.95)

        results = {
            "node 1":node_1.name,
            "GDI 1": np.round(node_1_gdi_mean, 2),
            "lower bound 1" : np.round(node_1_gdi_hpd[0],2),
            "upper bound 1" : np.round(node_1_gdi_hpd[1],2),
            
            "node 2":node_2.name,
            "GDI 2": np.round(node_2_gdi_mean,2),
            "lower bound 2" : np.round(node_2_gdi_hpd[0],2),
            "upper bound 2" : np.round(node_2_gdi_hpd[1],2),
            
            f"{cf_dict['mode']} accepted?": node_1.modified}

        return results

    # collect feedback
    feedback = []
    for pair in node_pairs_to_modify:
        feedback.append(node_pair_feedback(pair[0], pair[1], gdi_values, cf_dict))

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
        gdi_values: Dict[NodeName, NumericParam],
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
        node_pair_decision(pair[0], pair[1], gdi_values, cf_dict)
    
    print_decision_feedback(node_pairs_to_modify, tree, gdi_values, cf_dict)

    # output the current imap of accepted species
    result_imap = get_attribute_filtered_imap(tree, attribute='species')
    imapfile_write(result_imap, 'RESULT_IMAP.txt')

    # output the currently accepted newick tree
    result_tree = get_attribute_filtered_tree(tree, attribute='species')
    with open('RESULT_TREE.txt', 'w') as f:
        f.write(result_tree)

    return tree