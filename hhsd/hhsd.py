'''
MAIN HHSD SCRIPT
'''

from sys import argv
import sys

from .customtypehints import CfileParam, BppCfileParam, Cfile
from .module_HA import HA_iteration, check_contintue, set_starting_state
from .module_bpp import bppctl_init
from .module_tree import init_tree
from .module_cmdline import cmdline_init
from .module_msa_imap import imapfile_read
from .module_cf_ingest import ingest_cf
from .module_helper import output_directory


# main wrapper function implementing pipeline functions
def hhsd(
        cf_path: Cfile,
        cf_override
        ):

    # read control file
    cf:CfileParam = ingest_cf(cf_path, cf_override)

    # set up the output directory
    output_directory(cf['output_directory'])

    # intialise bpp control file 
    bpp_ctl:BppCfileParam = bppctl_init(cf)

    # read in essential data
    imap = imapfile_read(imap_filename=cf['Imapfile'], output_type="popind")
    newick = cf['guide_tree']

    # initailise tree
    tree = init_tree(newick, imap)

    # set up the starting proposal
    tree = set_starting_state(tree, cf['mode'])

    # run iterative algorithm
    while True:
        
        tree = HA_iteration(tree, bpp_ctl, cf)
        
        if check_contintue(tree, cf):
            continue
        else:
            sys.exit("Quitting hhsd")


"""
MAIN ENTRY POINT FOR RUNNING HHSD FROM THE TERMINAL
"""
def run():
    cf_path, cf_override = cmdline_init(argv)
    hhsd(cf_path, cf_override)