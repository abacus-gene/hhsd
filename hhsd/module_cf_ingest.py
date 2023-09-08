'''
FUNCTIONS REQUIRED TO INGEST AND CHECK THE PARAMETERS OF THE CONTROL FILE PROVIDED MY THE USER.
'''

import pandas as pd
import io
import sys
from collections import Counter

from .module_helper import read_filter_comments, read_format_curlybrackets, stripall, dict_merge, param_name_match
from .module_check_helper_cf import check_output_dir, check_msa_file, check_imap_file, check_newick, check_imap_msa_compat, check_imap_tree_compat, check_can_infer_theta, check_mode, check_GDI_threshold, check_migration
from .module_check_helper_bpp import check_seed, check_tauprior, check_thetaprior, check_finetune, check_sampfreq, check_nsample, check_burnin, check_locusrate, check_cleandata, check_threads, check_threads_msa_compat, check_nloci, check_nloci_msa_compat, check_threads_nloci_compat, check_migprior, check_phase

# dictionary of CF parameters that are currently supported
cf_param_dict =    {
    # location of output 
    "output_directory"      :None,
    # data sources and empirical parameters
    "seqfile"               :None,
    "Imapfile"              :None,
    "guide_tree"            :None,
    "phase"                 :None,
    #"mrate"                 :None,
    # parameters for the hierarchical method
    "mode"                  :None,
    "GDI_threshold"          :None,
    #"generation_threshold"  :None,
    #"decision_criteria"     :None,
    # parameters passed to BPP instances
    "seed"                  :None,
    "thetaprior"            :None,
    "tauprior"              :None,
    "finetune"              :None,
    "sampfreq"              :None,
    "nsample"               :None,                   
    "burnin"                :None,
    "threads"               :None,
    "nloci"                 :None,
    "locusrate"             :None,
    "cleandata"             :None,
    # migration related parameters
    "migprior"              :None,
    "migration"             :None,
}


## FUNCTIONS FOR READING AND PROCESSING THE CF INTO PYTHON COMPATIBLE OBJECTS

# read the control file into a dataframe
def read_cf_to_df(
        cf_name,
        ):
    
    # strip out all comments from control file
    try:
        lines = read_filter_comments(cf_name)
        lines_text = "\n".join(lines)
    except:
        sys.exit("Error: could not remove comments from control file.\nCheck formatting and refer to manual.")
    
    # format parameters bounded by '{}'
    try:
        lines_text = read_format_curlybrackets(lines_text)
    except:
        sys.exit("Error: could not parse control file parameters delimited by '{' and '}'.\nCheck formatting and refer to manual.")

    # read into dataframe
    try:
        buf = io.StringIO(lines_text)
        df = pd.read_csv(buf, sep = "=", header=None)
        
        # rename columns, convert every cell to text, and remove trailing and leading whitespaces
        df.rename(columns={0: "par", 1: "value"}, inplace=True)
        for col in df.columns:
            df[col] = df[col].map(stripall)
    except:
        sys.exit("Error: could not parse control file\nCheck formatting and refer to manual.")
    
    return df

# transform the dataframe into a dict 
def cf_df_to_dict(
        cf_df
        ):

    # make dataframe a dict, and transform all types to string
    cfdict = cf_df.set_index('par').T.to_dict("records")[0]
    cfdict = {str(param):str(cfdict[param]) for param in cfdict}

    # if a given parameter is present, but but has an empty value, return None for that parameter
    for param in cfdict:
        if len(cfdict[param]) == 0:
            cfdict[param] = None

    cfdict = dict_merge(cf_param_dict, cfdict)

    return cfdict

## FUNCTION FOR CHECKING PARAMETER NAMES
def verify_cf_parameter_names(
        cf_file,
        cf_override
        ):

    '''
    This function aims to ensure that the control file only contains calls
    to valid parameters of the pipeline, and does not have duplicate values. 
    '''

    cf_df = read_cf_to_df(cf_file)
    
    # check if all parameters are from the known correct list
    supplied_param = cf_df['par'].tolist()
    correct_param = list(cf_param_dict.keys())

    unmatched = set(supplied_param).difference(set(correct_param))
    if len(unmatched) > 0:
        error_msg = "Error: unknown parameters in control file:"
        for param in unmatched: error_msg += f"\n\t{param} {param_name_match(param, correct_param)}"
        sys.exit(error_msg)
    
    # check for duplicate values
    counts = dict(Counter(supplied_param))
    if any(value > 1 for value in counts.values()):
        error_msg = "Error: duplicate parameters in control file:\n\t"
        error_msg += str([parameter for parameter in counts if counts[parameter] > 1])[1:-1]
        error_msg += "\ndelete or comment out lines with duplicate parameters"
        sys.exit(error_msg)

    # return cf dict if all the tests were passed
    cf = cf_df_to_dict(cf_df)

    # TEMPORARY SOLUTION FOR OVERRIDE, FIX LATER TO ADD CHECKS
    if cf_override != None:
        cf = dict_merge(cf, cf_override)

    return cf



####
'''
Checks all values provided in a control file to ensure that the program will not crash.
If any of these checks do fail, the program will quit with an informative error message. 
'''
def cf_parameter_check(
        cf
        ):

    #  Checking parameters of the control file (functions explained and implemented in 'module_check_helper_cf')
    cf['output_directory'] = check_output_dir(cf['output_directory'])

        # check data is of correct type
    cf['seqfile']  = check_msa_file(cf['seqfile'])
    cf['Imapfile'] = check_imap_file(cf['Imapfile'])
    check_newick(cf['guide_tree'])
    
        # compatibility checking of data
    check_imap_msa_compat(cf['Imapfile'], cf['seqfile'])
    check_imap_tree_compat(cf['Imapfile'], cf['guide_tree'])
    check_phase(cf['phase'])
    check_can_infer_theta(cf['Imapfile'], cf['seqfile'], cf['phase'],)
    
        # check parameters of the hierarchical method
    check_mode(cf['mode'])
    cf['GDI_threshold'] = check_GDI_threshold(cf['GDI_threshold'], cf['mode'])

    # Checking parameters passed to BPP(functions explained and implemented in 'module_check_helper_bpp')
    check_seed(cf['seed'])
    check_tauprior(cf['tauprior'])
    check_thetaprior(cf['thetaprior'])
    check_finetune(cf['finetune'])
    check_sampfreq(cf['sampfreq'])
    check_nsample(cf['nsample'])
    check_burnin(cf['burnin'])
    check_locusrate(cf['locusrate'])
    check_cleandata(cf['cleandata'])
    
    check_threads(cf['threads'])
    check_threads_msa_compat(cf['threads'], cf['seqfile'])
    
    check_nloci(cf['nloci'])
    check_nloci_msa_compat(cf['nloci'], cf['seqfile'])
    check_threads_nloci_compat(cf['threads'], cf['nloci'])

    check_migprior(cf['migprior'])

    # Checking parameters related to migration,
    cf['migration'] = check_migration(cf["migration"], cf['migprior'], cf['guide_tree'])
    

    return cf


## FINAL WRAPPER FUNCTION IMPLEMENTING READING AND CHECKING
def ingest_cf(
        cf_file,
        cf_override
        ):

    # verify that the cf only contains known parameters, and no duplicates
    cf = verify_cf_parameter_names(cf_file, cf_override)
    
    # verify the specific values provided for the parameters will allow the program to run
    cf = cf_parameter_check(cf)

    return cf