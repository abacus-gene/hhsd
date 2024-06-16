'''
FUNCTIONS RELATED TO READING, MODIFYING, AND WRITING TREE DATA STRUCTURES
'''

import copy
from typing import Literal, Union, Tuple, List, Dict

from .customtypehints import NodeName, NewickTree, ImapPopInd, ImapIndPop, AlgoMode
from .module_ete3 import Tree, TreeNode


def tree_to_newick  (
        tree:           Tree
        ) ->            NewickTree:
    
    '''
    small wrapper function that returns an ete3 tree in a newick formatted string
    '''

    return tree.write(format=9)

## FUNCTIONS FOR INITIALIZING THE MAIN TREE OBJECT BASED ON THE DATA PROVIDED BY THE USER
def name_internal_nodes(
        tree:   Tree
        ) ->    Tree:

    '''
    Takes a tree where only leaf nodes are named, and names all internal nodes. \\
    This is accomplished by combining the names of the descendant nodes. \\
    This is done from leaf to root. e.g. leaf nodes 'A' & 'B' will have ancestor 'AB'
    '''

    for node in tree.traverse("postorder"):
        if len(node.name) == 0:
            newname = ""
            for count, descendants in enumerate(node.iter_descendants("levelorder")):
                if count == 0 or count == 1:
                    newname += descendants.name
            node.name = newname
    
    return tree

def add_individuals_to_tree(
        pop_tree:       Tree,
        imap_popind:    ImapPopInd,
        ) ->            Tree:
    
    '''
    Based on the imap, add individuals from each population as child nodes (with attribute "individual_node") to the leaf population nodes. 
    '''
    
    # iterate through each population
    for pop_name in imap_popind:
        # begin by finding the guide tree node with the name matching the population
        pop_node = list(filter(lambda n: n.name == pop_name and n.node_type == "population", pop_tree.traverse()))[0]
        
        # iterate through each individual in the population, and add as a child
        for ind_name in imap_popind[pop_name]: 
            pop_node.add_child(name = ind_name)
        
        # add child attribute to resulting nodes
        for child_node in pop_node.iter_descendants(): 
            child_node.add_features(node_type = "individual",)

    return pop_tree


def init_tree(
        tree_newick:    NewickTree,
        imap_popind:    ImapPopInd,
        ) ->            Tree:

    '''
    Initialize the Tree object, based on the 'guide_tree' parameter, 
    and add individuals from the 'Imapfile' onto the tree as child nodes of their parent populations.
    '''

    # ingest newick and turn into ete3 tree
    tree_out = Tree(newick = tree_newick)
    
    # name internal nodes based on leaf nodes
    tree_out = name_internal_nodes(tree_out)

    # set all nodes to be populations 
    for node in tree_out.traverse("postorder"): node.add_features(node_type = "population",)

    # add iteration attribute to root node, which is initally set to 0
    root = tree_out.get_tree_root(); root.add_features(iteration = 0)

    # add nodes corresponding to individuals to the tree
    tree_out = add_individuals_to_tree(tree_out, imap_popind)

    return tree_out

## FUNCTIONS FOR OUTPUTTING DIFFERENT TYPES OF DATA BASED ON THE CURRENT STATE OF THE TREE

def get_attribute_filtered_tree(
        tree:           Tree,
        attribute:      Literal['merge','split','population'],
        newick =        True, #returns a newick representation of the filtered tree by default, but can also be used to return the full tree object
        ) ->            Union[NewickTree, Tree]:
    
    '''
    This function is used to filter the full tree datastructure to only include certain nodes. Use cases include:
    1) Filtering for only accepted species, which can be done at the end of an iteration
    2) Filtering nodes such that the tree used to evaluate split proposals is output, this is done at the begenning of a split iteration
    3) Filtering for populations only. This is useful whenever the guide tree needs to be accessed
    '''

    copy_tree = copy.deepcopy(tree)

    # filter to currently accepted species. This is also utilized as the proposal for merge mode
    if attribute == "species" or attribute == "merge":
        nodes = copy_tree.search_nodes(species = True)
                         
        # base case
        if len(nodes) != 1:
            copy_tree.prune(nodes)
        # special exception case when only the root node is a species, 
        # as ete3 pruning would return the whole tree, which is not the intended behaviour
        else: 
            root = copy_tree.get_tree_root()
            descendants = list(root.iter_descendants("levelorder"))[0:2]
            copy_tree.remove_child(descendants[0])
            copy_tree.remove_child(descendants[1])
    
    # filter to currently accepted species + species resulting from a split proposal 
    elif attribute == "split":
        copy_tree.prune(list(filter(lambda n: n.proposal == "split" or n.species == True, copy_tree.search_nodes(node_type = "population"))))

    # filter to populations only
    elif attribute == "population":
        copy_tree.prune(copy_tree.search_nodes(node_type = "population"))

    if   newick == True:
        return tree_to_newick(copy_tree)
    elif newick == False:
        return copy_tree


def get_attribute_filtered_imap(
        tree:           Tree,
        attribute:      Literal['merge','split','species'],
        ) ->            ImapIndPop:
    
    '''
    For each individual in the input 'Imapfile', find the first population node relevant to the current problem. This can be:
    1) The population the individual would be assigned to under a merge or split proposal 
    2) The currently accepted species the individual is mapped to
    '''
    
    indpop_dict:ImapIndPop = {}
    
    # isolate individual nodes 
    for node in tree.search_nodes(node_type="individual"):
       
        parent_candidate = node.up
        
        # if looking for accepted species (which are also the merge candidates)
        if attribute == "species" or attribute == "merge":
            while parent_candidate.species == False:
                parent_candidate = parent_candidate.up
        
        # if looking for the nodes required to implement the split proposal
        elif attribute == "split":
            while parent_candidate.proposal != "split" and parent_candidate.species == False:
                parent_candidate = parent_candidate.up
        
        # add the parent for the specific individual node to the indpop dict
        indpop_dict[f"{node.name}"] = f"{parent_candidate.name}"

    return indpop_dict

def get_current_leaf_species(
        tree:           Tree,
        ) ->            List[NodeName]:

    '''
    Get the list of currently accepted species nodes that are also leaves.
    '''

    copy_tree = get_attribute_filtered_tree(tree, "species", newick=False)
    leaf_names = [leaf.name for leaf in copy_tree]

    return leaf_names

def get_iteration(
        tree:   Tree
        ) ->    int:
    
    '''
    get the iteration from the 'iteration' attribute located on the root node
    '''

    root = tree.get_tree_root()
    return root.iteration    

def get_all_populations(
        tree_newick:    NewickTree,
        ) ->            List[NodeName]:
    
    '''
    Get the list of all populations including internal nodes based on a newick tree.
    '''
        
    tree_out = Tree(newick = tree_newick)
    tree_out = name_internal_nodes(tree_out)
    population_names = [node.name for node in tree_out.traverse()]

    return population_names

def get_first_split_populations(
        tree_newick:    NewickTree,
        ) ->            Tuple[List[NodeName],List[NodeName]]:

    '''
    Get two lists for the guide tree leaf populations corresponding to the two first branches.
    '''

    tree = Tree(newick = tree_newick)
    root = tree.get_tree_root(); children = root.get_descendants()
    
    node_1 = children[0];                   node_2 = children[1]
    pop_1 = node_1.get_leaves();            pop_2 = node_2.get_leaves()
    pop_1 = [node.name for node in pop_1];  pop_2 = [node.name for node in pop_2]

    return pop_1, pop_2


def get_node_pairs_to_modify(
        tree:           Tree,
        algo_mode:      AlgoMode
        ) ->            List[Tuple[TreeNode, TreeNode]]:

    '''
    Get a list of node pairs to modify according to the merge or split proposal.
    '''

    # search depending on mode
    if   algo_mode == "merge":
        node_pairs_to_modify = tree.search_nodes(proposal='merge')
    elif algo_mode == "split":
        node_pairs_to_modify = tree.search_nodes(proposal='split')
    
    # split into pairs of nodes
    node_pairs_to_modify = [node_pairs_to_modify[x:x+2] for x in range(0, len(node_pairs_to_modify), 2)]

    return node_pairs_to_modify


## FUNCTIONS FOR ADDING NUMERICAL NODE ATTRIBUTES (E.G. TAU OR M) INFERRED USING BPP TO THE TREE OBJECT
def add_attribute_tau_theta(
        tree:           Tree,
        tau_values:     Dict[NodeName, float],
        theta_values:   Dict[NodeName, float],
        ) ->            Tree:

    '''
    Add tau and theta parameters inferred using A00 as node attributes.
    '''
    
    # add tau and theta attributes to all population nodes
    for node in tree.search_nodes(node_type="population"):
        # tau values are only available for non-leaf (ancestral) nodes.
        if node.name in list(tau_values.keys()):
            node.add_features(tau = float(tau_values[node.name]))
        else:
            node.add_features(tau = None)

        # theta values should be available for all nodes.
        if node.name in list(theta_values.keys()):
            node.add_features(theta = float(theta_values[node.name]))
        else:
            node.add_features(theta = None)

    return tree