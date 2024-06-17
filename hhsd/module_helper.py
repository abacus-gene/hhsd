'''
HELPER FUNCTIONS USED BY OTHER MODULES
'''

import re
import copy
import sys
import os
from pkg_resources import resource_filename
import platform
from pathlib import Path
from difflib import SequenceMatcher

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

# read a text file, and return a filtered version with all text after "*" removed
def read_filter_comments(
        input_file
        ) -> list[str]:

    lines = readlines(input_file)
    lines = [line.split("#")[0] for line in lines]
    lines = [line.split("*")[0] for line in lines]
    lines = remove_empty_rows(lines)

    return lines

# format any parameters that use '{' '}' to spread over multiple lines into a single line
def read_format_curlybrackets(
        text_lines: str
        ) -> str:
    
    # split into normal sections, and sections bounded by '{' '}'
    temp_sections = re.split(r'([{][^\{\}]*[}])', text_lines)
    temp_sections = [section for section in temp_sections if len(section) > 0]
    
    # in the sections bounded by '{' '}', remove any newline characters
    temp_sections = [re.sub('\n','', item) if (item[0] == "{" and item[-1] == "}") else item for item in temp_sections]
    
    return "".join(temp_sections)

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


# check if a single numeric parameter is supplied in the correct format and range
def check_numeric(
        value, 
        statement = "-10000000<x<1000000", 
        float_or_int = "f",
        ):
    
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
    get the correct platform-specific path to the bundled bpp executable, depending on the platform of the user
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
        exec_path = os.path.join(bpp_folder, 'macos', 'bpp')
    
    else:
        sys.exit(f"HHSD does not support the current operating system: {system}")

    return exec_path