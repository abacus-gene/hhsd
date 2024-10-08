'''
UNCTIONS REQUIRED FOR INTERACTING WITH THE PROGRAM FROM THE COMMAND LINE
'''

import os
from pathlib import Path
from pathlib import PurePath
import sys
import multiprocessing

from .customtypehints import Cfile, CfileParam
from .module_helper import stripall


# move to the directory where the pipeline is located.
def resolve_cf_file(
        filepath: Cfile
        ):

    '''
    Moving to the working directory is crucial, as filepaths in the control file are interpreted relative to
    the working directory. Output from the pipeline will also be located in the working directory.
    '''

    # check if the path provided exists
    try:
        cf_filepath = Path(filepath)
        cf_filepath = cf_filepath.resolve(strict=True)
    except:
        sys.exit(f"FilePathError: file path '{filepath}' of control file could not be resolved.")

    # check the directory where the cf is supposed to be located, and move in if possible
    try: 
        cf_directory = PurePath(cf_filepath).parent
        os.chdir(cf_directory)
    except:
        sys.exit(f"FilePathError: could not access requested folder of control file located at: '{filepath}'.")
    
    try:
        cf_filename = PurePath(cf_filepath).name
        return cf_filename
    except:
        sys.exit(f"MissingControlFileError: no control file present at file present at: '{filepath}'")


# used to separate command line arguments into categories
def group(
        seq, 
        sep
        ):
        
    g = []
    for el in seq:
        if el in sep:
            yield g
            g = []
        g.append(el)
    yield g

def categorise_arguments(
        argument_list
        ):
    
    # separate commands into categories
    argument_categories = list(group(argument_list, ['--cfile','--cfpor']))[1:]
    # get the string of the parameters in a non-empty category
    argument_categories = {cat[0]:" ".join(cat[1:]) for cat in argument_categories if len(cat) > 1 } 

    return argument_categories

def interpret_parameter_override(
        argument
        ) -> CfileParam:

    args_to_override = argument.split(", ")
    
    try:
        cf_ovveride_dict = {stripall(arg.split("=")[0]):stripall(arg.split("=")[1]) for arg in args_to_override}
    except:
        sys.exit(f"ParameterOverrideError: --cfpor argument '{argument}' used incorrect syntax")

    return cf_ovveride_dict

## FINAL WRAPPER FUNCTION
def cmdline_init(
        argument_list
        ):

    # check that the required arguments are provided
    arguments_dict = categorise_arguments(argument_list[1:])

    # load splash text if no arguments are provided
    if len(arguments_dict) == 0 and "--cfile" not in argument_list:
        sys.exit(f"hhsd version 1.0.0\n{multiprocessing.cpu_count()} cores available\nspecify control file for analysis with --cfile")

    # check that the control file is specified
    if "--cfile" not in arguments_dict:
        sys.exit("MissingControlFileError: please specify control file as '--cfile name_of_control_file'")
    cf_path = arguments_dict['--cfile']
    cf_path:Cfile = resolve_cf_file(cf_path) # move to the directory where the cf is located

    print("\n< Checking control file arguments... >\n")

        # check any parameter overrides
    if "--cfpor" in arguments_dict:
        cf_override_dict:CfileParam = interpret_parameter_override(arguments_dict['--cfpor'])
        print("The following control file parameters have been overridden via --cfpor:")
        for element in cf_override_dict: print(element, "=", cf_override_dict[element])
        print("\n")
    else:
        cf_override_dict = None

    return cf_path, cf_override_dict