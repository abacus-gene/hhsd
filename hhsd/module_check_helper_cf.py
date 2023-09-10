'''
FUNCTIONS THAT AIM TO CHECK FOR  MISSPECIFICATIONS OF DATA PROVIDED IN THE 
CONTROL FILE, AND PARAMETERS RELEVANT TO THE HM ALGORITHM
'''

import sys
import re
from pathlib import Path

from .module_ete3 import Tree
from .module_helper import check_file_exists, check_folder
from .module_msa_imap import alignfile_to_MSA, imapfile_read, count_seq_per_pop
from .module_tree import name_internal_nodes, get_all_populations
from .module_migration import read_specified_mig_pattern

# check if the output directory is available
def check_output_dir(
        output_directory
        ):

    if output_directory == None:
        sys.exit("cfile error: no 'output_directory' provided. specify a folder where the results of the analysis should be deposited")
    check_folder(output_directory)
    
    final_output_directory = Path(output_directory).resolve(strict=False)
    if str(final_output_directory) != output_directory:
        print(f"filepath for output direcectory inferred to be:\n\t{final_output_directory}")

    return final_output_directory

## FILE TYPE CHECKS
# check if an alignment file can be loaded in as a valid MSA object
def check_msa_file(
        seqfile
        ):

    if seqfile == None:
        sys.exit("cfile error: 'seqfile' not specified")
    
    check_file_exists(seqfile, 'seqfile')
    
    # try to load the alignment file to the internal MSA object
    try:
        align = alignfile_to_MSA(seqfile)
    except:
        sys.exit(f"cfile error: 'seqfile' {seqfile} is not a valid phylip MSA")

    # check that all sequence ids are formatted correctly
    for curr_locus in align:
        curr_id_list = [seq.id for seq in curr_locus]
        for id in curr_id_list:
            if not bool(re.fullmatch(r"^\S+\^\S+$|^\^\S+$", id)):
                sys.exit(f"Error: sequence names in 'seqfile' '{seqfile}' do not follow requred naming conventions. \nSequence names should be in the format seq_id^individual_id or ^individual_id.")

    final_seqfile = Path(seqfile).resolve(strict=True)

    if str(final_seqfile) != seqfile:
        print(f"filepath for seqfile inferred to be:\n\t{final_seqfile}")

    return final_seqfile

    
# check if the file supposted to be an imap is actually an Imap
#### FIX FIX NEEDS EXTRA WORK
def check_imap_file(
        imapfile
        ):

    if imapfile == None:
        sys.exit("cfile error: 'Imapfile' not specified")
    
    check_file_exists(imapfile, 'imapfile')

    # try to load the imap in both modes
    try:
        imap = imapfile_read(imapfile, "popind")
        imap = imapfile_read(imapfile, "indpop")
    except:
        sys.exit(f"cfile error: 'Imapfile' {imapfile} is not of the correct filetype")

    final_imapfile = Path(imapfile).resolve(strict=True)

    if str(final_imapfile) != imapfile:
        print(f"filepath for Imapfile inferred to be:\n\t{final_imapfile}\n")

    return final_imapfile

# check if a supplied tree is correctly, formatted, and contains no polytomies
def check_newick(
        tree
        ):

    # check if tree is supplied    
    if tree == None:
        sys.exit("cfile error: 'guide_tree' not supplied")
    
    # check if tree can be ingested as a newick tree
    else:
        try:
            t = Tree(tree)
        except:
            sys.exit("cfile error: 'guide_tree' could not be processed. Check formatting.")
    
    
    node_names = [node.name for node in t.traverse("postorder") if node.name != ""]

    # check for less than 2 nodes
    if len(node_names) < 2:
        sys.exit("Error: 'guide_tree' must have at least two nodes")

    # check for repeated node names
    if len(node_names) > len(set(node_names)):
        sys.exit("Error: 'guide_tree' has repeated node names")

    # check for conflicts arising from overlaps between leaf and node names
    t = name_internal_nodes(t)
    node_names = [node.name for node in t.traverse("postorder")]
    if len(node_names) > len(set(node_names)):
        sys.exit("Error: 'guide_tree' internal node names overlap with leaf node names")

    # check if the tree is binary, and return a special error if not
    for node in t.traverse():
        if len(node.get_children()) not in [0, 2]:
            sys.exit("Error: 'guide_tree' is not binary")

    return True


## DATASET COMPATIBILITY CHECKING
# check if an Imap file and a sequence alignment are mutually compatible
'''
This fucntion checks that all individual IDs in the alignment are mapped 
to a population in the IMAP. If this is the case, BPP can proceed successfullu
'''
def check_imap_msa_compat(
        imapfile, 
        seqfile
        ):

    # get the list of individual IDs in the Imap
    names_imap = set(list(imapfile_read(imapfile, "indpop").keys()))
    
    # get the list of individual IDs in the alignment
    alignment = alignfile_to_MSA(seqfile)
    names_align = set()
    for locus in alignment:
        for seq in locus:
            name = seq.id
            name = name.split("^")[-1]
            names_align.add(name)

    # check if the two sets of names are not identical
    if names_imap != names_align:
        
        error_msg = "cfile error: Imap and seqfile are incompatible\n"
        
        not_found_in_alignment = names_imap.difference(names_align)
        if len(not_found_in_alignment) > 0:
            error_msg += "\n\tthe following IDs are in the Imap, but not the seqfile:\n"
            error_msg += f"\t{str(not_found_in_alignment)[1:-1]}\n"
        
        not_found_in_imap = names_align.difference(names_imap)
        if len(not_found_in_imap) > 0:
            error_msg += "\n\tthe following IDs are in the seqfile, but not the Imap:\n"
            error_msg += f"\t{str(not_found_in_imap)[1:-1]}\n"
        
        sys.exit(error_msg)

# check if a Newick tree and an Imap file are mutually compatible
'''
This implies that the all of the tree leaf names in the newick string are also
found in the Imap file, and vice versa.
'''
def check_imap_tree_compat(
        imapfile, 
        tree
        ):

    # get the list of population names in the Imap
    pops_imap = set(list(imapfile_read(imapfile, "popind").keys()))
    
    # get the list of populations mentioned in the tree
    pops_tree = set()
    t = Tree(tree)
    for node in t.traverse("levelorder"):
        if node.name != "": pops_tree.add(node.name)
    
    # if not, provide detailed feedback about the missing populations
    if pops_imap != pops_tree:
        
        error_msg = "cfile error: Imap and guide tree are incompatible\n"
        
        not_found_in_alignment = pops_imap.difference(pops_tree)
        if len(not_found_in_alignment) > 0:
            error_msg += "\n\tthe following populations are in the Imap, but not the guide tree:\n"
            error_msg += f"\t{str(not_found_in_alignment)[1:-1]}\n"
        
        not_found_in_imap = pops_tree.difference(pops_imap)
        if len(not_found_in_imap) > 0:
            error_msg += "\n\tthe following populations are in the guide tree, but not the Imap:\n"
            error_msg += f"\t{str(not_found_in_imap)[1:-1]}\n"
        
        sys.exit(error_msg)

# check if a guide tree and Imap and alignment together are suitable for GDI calculations
'''
This means that each population in the starting imap has at least two sequences. Theta cannot be
esimated if there are less than two sequences, so this is required for gdi calculations
'''
def check_can_infer_theta(
        imapfile, 
        seqfile,
        phase,
        ):
    
    indpop_dict = imapfile_read(imapfile, "indpop")
    alignment = alignfile_to_MSA(seqfile)

    # count the max number of sequences per population
    seq_per_pop = count_seq_per_pop(indpop_dict, alignment)
    
    if phase == None or phase == "0":
        # if all sequences are unphased, then two sequences per population are required to estimate theta
        insufficient = [pop for pop in seq_per_pop if seq_per_pop[pop] < 2]

        if len(insufficient) > 0:
            error_msg = "cfile error: insufficient number of sequences in:\n"
            for pop in insufficient: error_msg += f"\t'{pop}' {seq_per_pop[pop]}\n"
            error_msg += "\nTheta cannot be estimated for unphased populations with 1 sequence.\nadd more sequences, remove species from the analysis, or specify phasing"
            sys.exit(error_msg)


## FUNCTIONS FOR SPECIFIC MCF PARAMETERS RELATED TO THE HM ALGORITHM

# check if the mode is correctly set to merge or split
def check_mode(
        mode
        ):

    if mode not in ['merge', 'split']:
        sys.exit("cfile error: please specify 'mode' as 'merge' or 'split'")

# check if the gdi threshold is correctly specified
def check_GDI_threshold(
        gdi_thresh,
        mode
        ):
        
    if gdi_thresh == None:
        sys.exit("cfile error: 'GDI_threshold' not specified. Please specify relation and cutoff value")
    else:
        if not bool(re.fullmatch("\A[<>]{1}[0]{1}[.]{1}[0-9]+\Z", gdi_thresh)):
            sys.exit(f"cfile error: 'GDI_threshold' incorrectly specified as '{gdi_thresh}'.\nPlease specify as '<0.x' or  '>0.x' (e.g. '<0.5' or '>0.8') depending on the mode")

        elif mode == "merge" and gdi_thresh[0] != "<":
            sys.exit(f"cfile error: in merge mode, the 'GDI_threshold' is an upper bound. specify as '<{gdi_thresh[1:]}' instead of '{gdi_thresh}'")
        elif mode == "split" and gdi_thresh[0] != ">":
            sys.exit(f"cfile error: in merge mode, the 'GDI_threshold' is a lower bound. specify as '>{gdi_thresh[1:]}' instead of '{gdi_thresh}'")

        return float(gdi_thresh[1:]) # return the parameter in the correct type


## MIGRATION SPECIFIC CHECKS

# in cases where migration is specified between certain populations, check that this is done correctly
def check_migration_newick_compatibility(
        mig_df, 
        tree_newick
        ):

    source_dest_df = mig_df.filter(items=['source', 'destination'])

    # check that there are no duplicates
    if any(mig_df.duplicated()):
        sys.exit("cfile error: same migration event specified more than once")

    # check if any of the migration events lack a source or destination
    missing_values_df = source_dest_df[source_dest_df.isnull().any(axis=1)]
    if len(missing_values_df["source"]) > 0:
        sys.exit(f"cfile error: the following specified migration events lack a source and/or destination:\n{missing_values_df.to_string()}")

    # check if all of the values correspond to known populations
    all_known_populations = get_all_populations(tree_newick)
    unknown_names_df = source_dest_df[source_dest_df[source_dest_df.isin(all_known_populations)].isnull().any(axis=1)]
    if len(unknown_names_df["source"]) > 0:
        sys.exit(f"cfile error: the following specified migration events have a source and/or destination population that is not found in the guide tree:\n{unknown_names_df.to_string()}")

    
    # check that no migration events occur between descendants 
    tree = Tree(newick = tree_newick)
    
    tree = name_internal_nodes(tree)

    for row in source_dest_df.itertuples(index = False):
        source_name = row[0]; dest_name =  row[1]
        source_node = tree.search_nodes(name = source_name)[0]; dest_node = tree.search_nodes(name = dest_name)[0]
        source_ancestors = source_node.iter_ancestors(); dest_ancestors = dest_node.iter_ancestors();
        if (source_node in dest_ancestors) or (dest_node in source_ancestors):
            sys.exit(f"cfile error: migration from '{source_name}' to '{dest_name}' not possible, as one is a descendant of the other.")


# check that the migration parameter has a valud value
def check_migration(
        migration,
        migprior,
        guide_tree
        ):
    
    # no migration
    if   migration == None:
        if migprior != None:
            sys.exit("cfile error: 'migprior' specified, but migration pattern was not.\nRemove migprior to analyse without migration, or specify migration patterns.")
        mig = None
    
    # specified migration
    elif migration[0] == "{" and migration[-1] == "}":
        if migprior == None:
            sys.exit("cfile error: migration pattern was specified, but 'migprior' not specified, but . Please specify a prior value for migration rates.")

        mig = read_specified_mig_pattern(migration)
        check_migration_newick_compatibility(mig, guide_tree)

    else:
        sys.exit("cfile error: migration parameter incorrectly specified. Refer to manual.")

    return mig



# # check if mutation rate is correctly specified
# def check_mrate(
#         mrate
#         ):
#     if mrate == None:
#         return None
#     else:
#         if not check_numeric(mrate, "0<x<1","f"):
#             sys.exit(f"cfile error: 'mrate' incorrectly specifed as '{mrate}'. use value between 0 and 1")
#         else:
#             return float(mrate)

# # check if generation threshold is correctly specified
# def check_generation_threshold(
#     gen_thresh,
#     mrate_status,
#     mode,
#     ):
    
#     if mrate_status == None and gen_thresh == None:
#         return None
#     elif mrate_status == True and gen_thresh == None:
#         sys.exit("Error: 'mrate' was specified, but 'generation_threshold' was not.")
#     else:
#         if mrate_status == None:
#             sys.exit("Error: 'generation_threshold' can only be set if 'mrate' is provided")
        
#         elif not bool(re.fullmatch("\A[<>]{1}[1-9]{1}[0-9]+\Z", gen_thresh)):
#             sys.exit(f"Error: 'generation_threshold' incorrectly specifed as '{gen_thresh}'.\nspecify as '<x' or '>x' (e.g. '<1000' or '>10000' depending on the mode)")
#         elif mode == "merge" and gen_thresh[0] != "<":
#             sys.exit(f"Error: in merge mode, the 'gen_threshold' is an upper bound. specify as '<{gen_thresh[1:]}' instead of '{gen_thresh}'")
#         elif mode == "split" and gen_thresh[0] != ">":
#             sys.exit(f"Error: in split mode, the 'gen_threshold' is a lower bound. specify as '>{gen_thresh[1:]}' instead of '{gen_thresh}'")
#         elif not 100 < int(gen_thresh[1:]) < 100000:
#             sys.exit("Error: 'gen_threshold' outside of permitted range. use value between 100 and 100000")

#         else:
#             return int(gen_thresh[1:]) # return the parameter in the correct type

# # check that the decision criteria is correctly set, and compatible with the other parameters
# def check_criteria(
#     criteria,
#     mrate_status,
#     gen_thresh_status
#     ):

#     # check if the decision criteria is in the accepted list
#     if criteria not in ['default', 'gdi_and_generations']:
#         sys.exit("Error: specify 'decision_criteria' as one of: 'any','any_two','all','age_&_one_gdi'")
    
#     # 'age_&_one_gdi' can only be used as a criteria, if both a mutation rate, and a generation threshold are specified
#     elif criteria == 'gdi_and_generations' and (mrate_status == None or gen_thresh_status == None):
#         sys.exit("Error: 'decision_criteria' is 'age_&_one_gdi' but 'mrate' and/or 'generation_threshold are not specifed.\nspecift 'mrate' and 'generation_threshold' or change criteria")