'''
FUNCTIONS REQUIRED TO INTERACT WITH BPP
'''

import subprocess
import re
import random
import sys
import time
import pandas as pd

from .customtypehints import CfileParam, BppCfileParam, BppCfile
from .module_helper import dict_merge, get_bundled_bpp_path, format_time
from .module_msa_imap import auto_prior, auto_nloci

# Custom styler function to format floats and tuples of floats with six decimal places
def format_float(value):
    if isinstance(value, float):
        return format(value, '.6f')
    elif pd.isna(value):
        return ' '  # Replace Pandas NaN with empty string
    else:
        return value


# contains the list of parameters that need to be present in a BPP control file
default_BPP_cfile_dict:BppCfileParam = {
    'seed':                 None,
    'seqfile':              None, 
    'Imapfile':             None, 
    'outfile':              "proposal_bpp_out.txt", 
    'mcmcfile':             "proposal_bpp_mcmc.txt",
    'speciesdelimitation':  '0',
    'speciestree' :         '0',
    'species&tree':         None, 
    'popsizes':             None, 
    'newick':               None,
    'phase':                '0',
    'usedata':              '1', 
    'nloci':                None,
    'locusrate':            '0',
    'cleandata':            '0', 
    'thetaprior':           None, 
    'tauprior':             None, 
    'finetune':             '1: .01 .0001 .005 .0005 .2 .01 .01 .01', 
    'print':                '1 0 0 0', 
    'burnin':               None, 
    'sampfreq':             '1', 
    'nsample':              None, 
    'threads':              None,
    'migprior':             None,
}

def bppctl_init(
        cf_param:   CfileParam
        ) ->        BppCfileParam:

    '''
    Runs at the start of the pipeline. Sets up BPP parameters that will be shared throughout the iterations.\\
    1) collects parameters of BPP that can be set through the control file (eg 'seed', 'threads'...)
    2) generates certain parameters that are needed for BPP to run, but not supplied.
    - seed is generated as a random int
    - nloci is the number of loci present in the MSA
    - tau and theta priors are inferred from the MSA and the delimitaiton, using the method of Bruce Rannala
    '''

    # ingest parameters from the mcf
    bpp_cdict:BppCfileParam = dict_merge(default_BPP_cfile_dict, cf_param)

    # generate values for parameters if not supplied
        # seed
    if bpp_cdict['seed'] == None:
        bpp_cdict['seed'] = random.randint(1, 100000)

        # nloci
    if bpp_cdict['nloci'] == None:
        bpp_cdict['nloci'] = auto_nloci(bpp_cdict['seqfile'])
        
        # tau and theta priors
    if bpp_cdict['tauprior'] == None or bpp_cdict['thetaprior'] == None:
        priors = auto_prior(imapfile=bpp_cdict['Imapfile'], seqfile=bpp_cdict['seqfile'], tree_newick=cf_param['guide_tree'])
        if bpp_cdict['tauprior']   == None:
            bpp_cdict['tauprior']   = priors['tauprior']
        if bpp_cdict['thetaprior'] == None:
            bpp_cdict['thetaprior'] = priors['thetaprior']


    return bpp_cdict

def bppcfile_write(
        bpp_param:      BppCfileParam, 
        ctl_file_name:  str
        ) ->            None: # writes bpp control file to disk
    
    '''
    Write dict representing the BPP ".ctl" file to a text file.
    '''
    # remove parameters without values
    bpp_param = {item:bpp_param[item] for item in bpp_param if bpp_param[item] != None}

    # convert to pandas dataframe
    df = pd.DataFrame(list(bpp_param.items()))
    
    # write dataframe to disk
    df.to_csv(ctl_file_name, sep = "=", header = False, index = False)
    with open(ctl_file_name, 'r+') as f:
        text = f.read()
        # substitute out these parameters, as they do not actually exist in bpp
        text = re.sub('popsizes=', '               ', text)
        text = re.sub('newick=', '               ', text)
        f.seek(0)
        f.truncate()
        f.write(text)


# run BPP with a given control file, and capture the stdout results
def run_BPP_A00(
        control_file:   BppCfile,  
        ) ->            None: # handles the bpp subprocess

    '''
    Handles the starting and stopping of the C program BPP, which is used to infer MSC parameters. 
    '''

    #subprocess.run(["bpp", "--cfile", control_file, "theta-move", "slide"]) # this line will make bpp fully verbose, only intended for debug purposes

    
    bpp_completed = False
    
    while not bpp_completed:
        start_time = time.time()
        progress = ' '

        # runs BPP in a dedicated subprocess
        process = subprocess.Popen(
            f"{get_bundled_bpp_path()} --cfile {control_file}", 
            shell = True, 
            bufsize = 1,
            stdout = subprocess.PIPE, 
            stderr = subprocess.STDOUT,
            encoding = 'utf-8', 
            errors = 'replace' 
            ) 
    
        # monitors the output of the process
        while True:
            time.sleep(0.01)
            realtime_output = process.stdout.readline()
            
            # check that numeric scaling is needed, and activate if it is.
            if "[ERROR] log-L for locus" in realtime_output:
                print('restarting BPP with numerical scaling', end = '\r')
                file1 = open(control_file, "a")  # append mode
                file1.write("scaling=1")
                file1.close()
                break
            
            # check if bpp errored with a given error message
            try:
                if "[ERROR]" in realtime_output:
                    sys.exit(f"Error: BPP failed. Check control file independently using the 'bpp' command")
            except:
                pass
        
            # check if bpp gave a segfault
            try:
                if "core dumped" in realtime_output:
                    sys.exit("Error: BPP failed with segfault. Check control file independently using the 'bpp' command")
            except:
                pass
                
            # proide feedback about completeness and time elapsed
            try:
                percentage = realtime_output.split()[0]
                if "%" in percentage and "." not in percentage: 
                    progress = percentage
            except:
                pass
            
            # print current state of progress
            print(f"BPP progress: {progress}      {format_time(time.time() - start_time)}                      ", end = "\r")
        
            # exit if process is stopped
            if realtime_output == '' and process.poll() is not None:
                print("                                                           ", end = "\r")
                return None
