
'''
THIS MODULE CONTAINS FUNCTIONS FOR VISUALIZING PROGRAM RESULTS
'''

import copy
import warnings
warnings.simplefilter('ignore')
import os

with warnings.catch_warnings():
    os.environ['QT_QPA_PLATFORM']='offscreen' # set qtl to nonwindowed mode, this way the pipeline should work through the command line

# EXTERNAL LIBRARY DEPENDENCIES
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree
    from ete3 import TreeStyle
    from ete3 import NodeStyle
    from ete3 import TextFace

import distinctipy

from module_tree import get_current_leaf_species

# make a tree with available distance values ultrametric, and scale the size to be easy to display
def make_ultrametric(tree):

    multiplier = (1/(tree.get_farthest_node()[1]))
    for node in tree.iter_descendants("postorder"):
        node.dist = node.dist*multiplier*3

    for node in tree.traverse("postorder"):
        descendants = list(node.iter_descendants("levelorder"))
        if len(descendants) > 2:
            node_1 = descendants[0]
            maxdist_1 = node_1.get_farthest_leaf()[1]
            node_2 = descendants[1]
            maxdist_2 = node_2.get_farthest_leaf()[1]
            if maxdist_1 > maxdist_2:
                node_2.dist = node_2.dist + (maxdist_1-maxdist_2)
            elif maxdist_2 > maxdist_1:
                node_1.dist = node_1.dist + (maxdist_2-maxdist_1)

    return tree


# graphically visualize the placement of individuals into species
def visualize_imap(tree, image_name = "imap.pdf"):
    
    tree = copy.deepcopy(tree)
    
    # get name of currently accepted populations
    curr_species = get_current_leaf_species(tree)

    # generate background colors corresponding to each currently accepted species
    colors = distinctipy.get_colors(len(curr_species), pastel_factor=0.5)
    colors = [distinctipy.get_hex(color) for color in colors]
    color_dict = {pop:colors[i] for i, pop in enumerate(curr_species)}

    # delete nodes corresponding to individuals
    for node in tree.traverse(): 
        if node.node_type == "individual":
            node.detach()

    # delete intermediate (non leaf) population nodes that are currently not accepted as species
    for node in tree.search_nodes(node_type="population"):
        if not node.is_leaf():
            if node.species == False:
                node.delete()

    # set up branch lengths according to available tau estimates
    if len(tree.search_nodes(species = True)) == 1: # this is only the case if the first split was rejected
        node = tree.get_tree_root()
        node.dist = node.tau
    else:
        for node in tree.iter_descendants():
            if node.species:
                anc = node.up
                node.dist = anc.tau
            else:
                node.dist = 0

    # set node style for ancestor nodes
    anc_style = NodeStyle()
    anc_style["size"] = 0
    anc_style["vt_line_width"] = 3
    anc_style["hz_line_width"] = 3
    anc_style["vt_line_color"] = "black"
    anc_style["hz_line_color"] = "black"

    for node in tree.traverse():
        if node.species == True and node.name not in curr_species:
            node.set_style(anc_style)
    
    # set node styles for accepted species and their descendants
    for node in tree.traverse():
        # if a node is a currently accepted species-
        if node.name in curr_species:
            descendants = list(node.iter_descendants())

            # add background color box to parent node, which unifies all individuals in a population into a single color
            sp_style = NodeStyle()
            sp_style["size"] = 0
            sp_style["vt_line_width"] = 0
            sp_style["hz_line_width"] = 3
            sp_style["vt_line_color"] = color_dict[node.name]
            sp_style["hz_line_color"] = "black"
            sp_style["bgcolor"] = color_dict[node.name]
            node.set_style(sp_style)
            
            # add the label of the species in the upper right hand corner of the color box
            label_node = list(node.traverse("postorder"))[0]
            face = TextFace(f"  {node.name} ", fsize=10, fstyle="bold")
            face.margin_left = 10
            face.hz_align = 2
            label_node.add_face(face, column=12, position = "aligned")
            #edited_parents.append(parent.name)

            # reduce the branch lengths of all descendant nodes to 0, and set their style such that they are invisible
            if len(descendants) > 0:
                for desc in descendants:
                    desc_style = NodeStyle()
                    desc_style["size"] = 0
                    desc_style["hz_line_width"] = 1
                    #desc_style["hz_line_color"] = color_dict[node.name]
                    desc.set_style(desc_style)
    
    # convert to ultrametric
    tree = make_ultrametric(tree)

    # set custom treestyle
    ts = TreeStyle()
    ts.branch_vertical_margin = 5
    ts.show_scale = False
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.scale = 100

    # final render
    with warnings.catch_warnings():
        tree.render(image_name, tree_style=ts)

# graphically visualize the placement of individuals into species
def visualize_decision(tree, image_name = "decision.pdf"):
    
    tree = copy.deepcopy(tree)
    
    # get name of currently accepted populations
    curr_species = get_current_leaf_species(tree)

    # delete nodes corresponding to individuals
    for node in tree.traverse(): 
        if node.node_type == "individual":
            node.detach()

    # delete intermediate any node that are not species or where not evaluated for species status
    for node in tree.traverse():
        if node.species == False and node.proposal == None:
            node.detach()

    # set up branch lengths according to available tau estimates
    for node in tree.iter_descendants():
        anc = node.up
        node.dist = anc.tau

    # set node style for ancestor nodes
    anc_style = NodeStyle()
    anc_style["size"] = 0
    anc_style["vt_line_width"] = 2
    anc_style["hz_line_width"] = 2
    anc_style["vt_line_color"] = "black"
    anc_style["hz_line_color"] = "black"

    for node in tree.traverse():
        node.set_style(anc_style)
    
    # add gdi and theta values to the nodes
    for node in tree.traverse():
        if node.proposal != None:
            face = TextFace(f"gdi: {node.gdi}\ntheta: {node.theta}", fsize=8)
            face.margin_left = 10
            face.hz_align = 2
            node.add_face(face, column=12, position = "aligned")

    # set up style for modified nodes
    mod_style = NodeStyle()
    mod_style["size"] = 0
    mod_style["vt_line_width"] = 4
    mod_style["hz_line_width"] = 4
    mod_style["vt_line_color"] = "green"
    mod_style["hz_line_color"] = "green"
    for node in tree.search_nodes(modified = True):
        node.set_style(mod_style)

    # set up style for unmodified nodes
    umod_style = NodeStyle()
    umod_style["size"] = 0
    umod_style["vt_line_width"] = 4
    umod_style["hz_line_width"] = 4
    umod_style["vt_line_color"] = "grey"
    umod_style["hz_line_color"] = "grey"
    for node in tree.search_nodes(modified = False):
        node.set_style(umod_style)
    
    # convert to ultrametric
    tree = make_ultrametric(tree)

    # set custom treestyle
    ts = TreeStyle()
    ts.branch_vertical_margin = 5
    ts.show_scale = False
    ts.margin_bottom = 10
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.scale = 100

    # final render
    with warnings.catch_warnings():
        tree.render(image_name, tree_style=ts)