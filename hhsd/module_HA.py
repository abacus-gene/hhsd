'''
EXECUTE A GIVEN ITERATION OF THE ITERATIVE DELIMITATION ALGORITHM
'''

import copy
import os

from .customtypehints import AlgoMode, CfileParam, BppCfileParam, MigrationPattern
from .module_ete3 import Tree
from .module_msa_imap import auto_pop_param, imapfile_write
from .module_helper import dict_merge
from .module_tree import get_attribute_filtered_imap, get_attribute_filtered_tree, get_current_leaf_species, add_attribute_tau_theta
from .module_bpp import bppcfile_write, run_BPP_A00, extract_param_estimate_from_outfile, extract_tau_theta_values, extract_mig_param_to_df
from .module_gdi_decision import tree_modify_delimitation, calculate_gdi
from .module_migration import append_migrate_rows


## MODIFICATION PROPOSAL RELATED FUNCTIONS


def set_starting_state(
        tree:   Tree,
        mode:   AlgoMode,
        ) ->    Tree:
    
    '''
    Runs before the iterations of the HM algorithm, and sets the starting species delimitation.
    - In merge mode, all populations in the guide tree start as species.
    - In split mode, only the root node is a species, as all populations are merged into a single species.

    Species status is marked with a node attribute "species".
    '''

    print('\n< Analysis started >\n')
    
    # in "merge" mode, all starting populations are initally accepted as species
    if mode == "merge":
        for node in tree.search_nodes(node_type="population"): node.add_features(species = True,)
    
    # in "split" mode, only the root node is initally accepted as species
    if mode == "split":
        for node in tree.search_nodes(node_type="population"): node.add_features(species = False,) 
        root_node = tree.get_tree_root()
        root_node.species = True

    # print feedback about starting state
    print(f"> Starting state of {mode} analysis\n")
    print(f"Number of species in starting delimitation:  {len(get_current_leaf_species(tree))}")
    print(str(get_current_leaf_species(tree))[1:-1])
    print(get_attribute_filtered_tree(tree, "species", newick=False))

    return tree


def set_tree_proposal_attributes(
        tree: Tree,
        mode: AlgoMode
        ) ->  Tree:
    
    '''
    Propose modifications to the topology by identifying the nodes whose species status will be assesed.
    - in merge mode, it indentifies leaf species node pairs that can be merged
    - in split mode, it identifies the descendant pairs of currently accepted leaf species nodes.

    The requirement for the species status of a node to be assesed is marked with the "proposal" attribute.
    '''

    tree = copy.deepcopy(tree)

    # start by marking all population nodes as non-leaf, and reset the status of current and accepted proposals
    for node in tree.search_nodes(node_type="population"): 
        node.add_features(leaf = False, proposal = None, modified = None)
    
    # mark leaf nodes in current delimitation
    for node in tree.search_nodes(species=True):
        # start my assuming a species node is a leaf
        leaf = True
        # check if it has descendants that are species, in which case it is not a leaf node
        descendants = node.iter_descendants()
        for descendant in descendants:
            if hasattr(descendant, "species"):
                if descendant.species == True:
                    leaf = False
        
        if leaf: node.leaf = True
    

    if   mode == 'merge':
    # mark merge candidates
        leaf_nodes = tree.search_nodes(leaf=True)
        for node in leaf_nodes:
            # check that the current node is not the root, as the root does not have ancestors
            if node != tree.get_tree_root():
                # get the current node and sister node as a list
                ancestor = node.up
                descendants = list(ancestor.iter_descendants("levelorder"))[0:2]
                # if both the current node and its sister are leaves, they are merge candidates
                if len(list(set(descendants) & set(leaf_nodes))) == 2:
                    for descendant in descendants:
                        descendant.proposal = "merge"
    

    elif mode == 'split':
    # mark split candidates
        leaf_nodes = tree.search_nodes(leaf=True)
        population_nodes = tree.search_nodes(node_type="population")
        for node in leaf_nodes:
            # get the two first degree children of the leaf node
            descendants = list(node.iter_descendants("levelorder"))[0:2]
            # if both children are populations (= not individuals), the children are split candidates
            if len(list(set(descendants) & set(population_nodes))) == 2:
                for descendant in descendants:
                    descendant.proposal = "split"


    return tree



def proposal_setup_files(
        tree:       Tree,
        bpp_cdict:  BppCfileParam,
        mode:       AlgoMode,
        migration:  MigrationPattern,
        ) ->        None: # writes files to disk
    
    '''
    Write the imap file and bpp control file needed to evaluate a given proposal
    '''

    # set up and write imap needed to evaluate proposal
    proposed_imap = get_attribute_filtered_imap(tree, mode)
    imapfile_write(proposed_imap, "proposed_imap.txt")
    
    # write control file needed to evaluate proposal
        # get proposed population parameters and topology
    prop_param = auto_pop_param(proposed_imap, bpp_cdict['seqfile'], bpp_cdict['phase'])
    prop_param['newick']    = get_attribute_filtered_tree(tree, mode)
    prop_param['Imapfile']  = "proposed_imap.txt"
    
    bpp_cdict = dict_merge(bpp_cdict, prop_param)
    bppcfile_write(bpp_cdict, "proposed_ctl.ctl")
    
    # if migration patterns are specified, append migration parameters to the control file
    if str(type(migration)) != "<class 'NoneType'>":
        append_migrate_rows(tree, migration, "proposed_ctl.ctl")


def HA_iteration(
        tree:       Tree, 
        bpp_cdict:  BppCfileParam, 
        cf_dict:    CfileParam
        ) ->        Tree:
    
    '''
    Important function implementing each iteration of the Hierarchical merge/split algorithm
    '''

    # increment iteration count
    root = tree.get_tree_root(); root.iteration = (root.iteration + 1)
    print(f"\n*** Iteration {root.iteration} ***")

    # create folder for iteration, and move in
    iter_dir_name = f"Iteration_{root.iteration}"
    os.mkdir(iter_dir_name)
    os.chdir(iter_dir_name)

    # inititate proposal by setting node attributes
    tree = set_tree_proposal_attributes(tree, cf_dict["mode"])

    # create bpp control file and imap file corresponding to proposal
    proposal_setup_files(tree, bpp_cdict, cf_dict["mode"], cf_dict["migration"])
    
    # run BPP and capture the output
    run_BPP_A00("proposed_ctl.ctl")
    estimated_param  = extract_param_estimate_from_outfile(BPP_outfile="proposal_bpp_out.txt")
    tau_values, theta_values = extract_tau_theta_values(estimated_param)
    migration_df     = extract_mig_param_to_df(estimated_param)

    # append the inferred numeric parameters to the tree
    tree = add_attribute_tau_theta(tree, tau_values, theta_values)

    # get gdi via calculations or simulations, and append results to the tree
    tree = calculate_gdi(tree, cf_dict['mode'], migration_df)

    # make decision based on results
    tree = tree_modify_delimitation(tree, cf_dict)

    # move back into working directory
    os.chdir("..")

    return tree    


# should the program move to the next HM iteration?
def check_contintue(
        tree:       Tree,
        cf_dict:    CfileParam
        ) ->        bool:

    '''
    Runs after each iteration of the hierarchical method, and checks if further iterations can be run, 
    or the algorithm has converged.
    '''

    next_iteration =  True

    # check if topology has been reduced to root node
    if   cf_dict['mode'] == 'merge':
        if len(get_attribute_filtered_tree(tree, "species", newick=False)) == 1:
            print("\nAll populations merged into single species. Final delimitation reached")
            next_iteration = False

    # check if topology has been expanded to the guide tree
    elif cf_dict['mode'] == 'split':
        if len(get_attribute_filtered_tree(tree, "population", newick=False)) == len(get_attribute_filtered_tree(tree, "species", newick=False)):
            print("\nAll populations in guide tree are species. Final delimitation reached")
            next_iteration = False
    
    # check if all propoosed modifications were rejected.
    if len(tree.search_nodes(modified=True)) == 0:
            print("\nAll modifications rejected. Final delimitation reached.")
            next_iteration = False

    return next_iteration