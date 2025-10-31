'''
HELPER FUNCTIONS USED BY OTHER MODULES
'''

import re
import copy
import sys
import os
import platform
from pathlib import Path
from difflib import SequenceMatcher
import subprocess

## CORE HELPER FUNCTIONS

# return a flattened list from a list of lists
def flatten (
        t: list[list]
        ) -> list:

    return [item for sublist in t for item in sublist]

# strip more than 1 trailing or leading whitespace
def stripall(
        input_string: str
        ) -> str:

    result = input_string
    try:
        while result[0] == " " or result[-1] == " ":
            result = result.strip()
    except:
        pass

    return result

# reads a text file into an array of rows
def readlines(
        file_name
        ) -> list[str]:

    fileObj = open(file_name, "r")
    lines = fileObj.read().splitlines() 
    fileObj.close()
    
    return lines 

# strip out any rows with only whitespace characters
def remove_empty_rows(
        input_rows: list[str]
        ) -> list[str]:

    output = [row for row in input_rows if re.search("\S+", row)]

    return output


# read values from one dictionary into the other if the value in the new dict is not none
def dict_merge(
        dict_1: dict,
        dict_2: dict
        ) -> dict:

    dict_out = copy.deepcopy(dict_1)

    # overwrite any valiues in the old dict where the new dict is not "?"
    for key in dict_out:
        if key in dict_2 and dict_2[key] != None:
            dict_out[key] = dict_2[key]

    return dict_out

# check if a supplied filename actually points to an existing file in the current directory
def check_file_exists(
        path,
        type_of_file = "file"
        ):

    # try to interpret the parameter as a path
    try:
        filepath = Path(path)
        filepath = filepath.resolve(strict=False) # will attempt to resolve the path, even if the file does not exist
    except:
        sys.exit(f"FilePathError: file path of {type_of_file} ('{path}') could not be resolved.")
    
    # check if location is a file
    if not filepath.is_file():
        sys.exit(f"FileNotFoundError: {type_of_file} location parameter '{path}' implies that \n'{filepath}' \nexists, but no such location was found in the file system.\nSpecify file paths relative to the directory of the control file, or use absolute paths.")    

# check if a specified folder does not exist, or is empty
def check_folder(
        path
        ):

    # try to interpret the parameter as a path
    try:
        filepath = Path(path)
        filepath = filepath.resolve(strict=False) # will attempt to resolve the path, even if the file does not exist
    except:
        sys.exit(f"FilePathError: file path of suggested 'output_directory' ('{path}') could not be resolved.")
    
    # check if location is a file
    if filepath.exists():
        if len(os.listdir(filepath)) != 0:
            sys.exit(f"ExistingFilesError: output directory '{path}' is non-empty.\ndelete files from '{filepath}'\nor provide the name of a new 'output_directory'")

def format_time(seconds):
    '''
    format time from seconds to H:M:S
    '''
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    return f'{round(hours)} h {round(minutes)} m {round(seconds)} s'



def check_numeric(
        value, 
        statement = "-10000000<x<1000000", 
        float_or_int = "f",
        ):
    """
    Check if a single numeric parameter is supplied in the correct format and range
    -value is a string corresponding to the number, e.g. "0.6", "8", or an incorrect input such as "number"
    -statement can be used to restrict the range of allowed values to some interval. 
    -float_or_int specifies if the number should be a float or int
    To be accepted, the input string has to be:
    1) a valid number (not other text)
    2) in the desired range (e.g. a probability cannot be less than 0 or greater than 1)
    3) the correct datatype (e.g. the number of MCMC samples cannot be 3.5)
    """
    
    # this module uses the literal eval function, which should only be used very very cautiously
    try: 
        # check the input conforms to the standard
        if not re.fullmatch("[-<>=x0-9.]+", statement):
            sys.exit("WARNING WARNING INCORRECT STATEMENT USED IN LITERAL EVAL")
    except:
        sys.exit("WARNING WARNING INCORRECT STATEMENT USED IN LITERAL EVAL")

    # actual checking of number
    try:
        # return true if number is of the right type, and in the correct range
        x = float(value)
        if float_or_int == "i":
            if bool(re.fullmatch("[0-9]+", value)) and x.is_integer() and eval(statement):
                return True
        elif float_or_int == "f":
            if bool(re.fullmatch("[0-9.]+", value)) and eval(statement):
                return True
    # return false otherwise
    except:
        return False

# get the closest matching parameter name from a list of valid values
def closest_param_match(
        paramname: str, 
        paramlist: list[str]
        ):

    match = {f"{element}":SequenceMatcher(None, paramname.lower(), element).ratio() for element in paramlist}
    match = [parameter for parameter in match if match[parameter] > 0.6 and match[parameter] == max(match.values())]
    if len(match) == 0:
        return ""
    else:
        return f"(Did you mean {str(match)[1:-1]}?)"


# handle the creation of the working directory if it is not already present
def output_directory(
        output_directory: Path
        ):

    if not output_directory.exists():
        output_directory.mkdir(parents = True)
    os.chdir(output_directory)


def get_bundled_bpp_path():
    '''
    get the correct OS-specific path to the bundled bpp executable, depending on the platform of the user
    '''
    
    # Determine the user's operating system
    system = platform.system()
    
    # Get the absolute path to the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define the absolute path from the script folder to the bpp folder
    bpp_folder = os.path.join(script_dir, 'bpp')

    if system   == 'Windows':
        exec_path = os.path.join(bpp_folder, 'windows', 'bpp.exe')
    
    elif system == 'Linux':
        exec_path = os.path.join(bpp_folder, 'linux', 'bpp')
   
    elif system == 'Darwin':  # macOS
        # Detect architecture: 'arm64' for Apple Silicon, 'x86_64' for Intel
        arch = platform.machine().lower()
        if arch == 'arm64':
            exec_path = os.path.join(bpp_folder, 'macos_arm', 'bpp')
        else:
            sys.exit(f"HHSD does not support MacOS computers with Intel processors. Please use a Linux or Windows machine, or a Mac with Apple Silicon.")
    
    else:
        sys.exit(f"HHSD does not support the current operating system: {system}")

    return exec_path

def check_bpp_executable():
    '''
    check that the bundled bpp executable is present and has executable permissions
    '''

    bpp_path = get_bundled_bpp_path()

    # check if the bpp executable exists
    if not os.path.isfile(bpp_path):
        sys.exit(f"BppExecutableError: the bundled bpp executable was not found at expected location: '{bpp_path}'")

    # check if the bpp executable has execute permissions
    if not os.access(bpp_path, os.X_OK):
        sys.exit(f"BppExecutableError: the bundled bpp executable at '{bpp_path}' does not have execute permissions. Please adjust the file permissions to allow execution.")

    # try to run the bpp command to check that it works
    try:
        result = subprocess.run([bpp_path, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            sys.exit(f"BppExecutableError: the bundled bpp executable at '{bpp_path}' could not be executed. Please ensure that the file is not corrupted and has the correct permissions.")
        if result.returncode == 0:
            return f"BPP executable at '{bpp_path}' is present and functional."
    except Exception as e:
        sys.exit(f"BppExecutableError: an error occurred while trying to execute the bundled bpp executable at '{bpp_path}': {e}")

    
