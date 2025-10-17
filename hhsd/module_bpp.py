'''
FUNCTIONS REQUIRED TO INTERACT WITH BPP
'''

import subprocess
import re
import random
import sys
import pandas as pd

from .customtypehints import CfileParam, BppCfileParam, BppCfile
from .module_helper import dict_merge, get_bundled_bpp_path
from .module_msa_imap import auto_prior, auto_nloci
from .module_tree import add_inner_node_names_to_newick

# contains the list of parameters that need to be present in a BPP control file
default_BPP_cfile_dict:BppCfileParam = {
    'seed':                 None,
    'seqfile':              None, 
    'Imapfile':             None, 
    'jobname':              'hhsd_job',
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
    'finetune':             '1', 
    'print':                '1 0 0 0', 
    'burnin':               None, 
    'sampfreq':             '1', 
    'nsample':              None, 
    'threads':              None,
    'wprior':               None,
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
        priors = auto_prior(
            imapfile=bpp_cdict['Imapfile'], 
            seqfile=bpp_cdict['seqfile'], 
            tree_newick=cf_param['guide_tree'],
            tau_prior=bpp_cdict['tauprior'], 
            theta_prior=bpp_cdict['thetaprior']
            )
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

    # add in internal node names to the newick string
    bpp_param['newick'] = add_inner_node_names_to_newick(bpp_param['newick'])

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
        ) ->            None: # handles the bpp subprocess, which outputs a file

    '''
    Handles the starting and stopping of the C program BPP, which is used to infer MSC parameters. 
    '''

    bpp_completed = False
    while not bpp_completed:
        # flag activated when numeric scaling is turned on
        restart = False

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
            output_line = process.stdout.readline()            
            
            # if the output is non-empty
            if output_line:
                # check that numeric scaling is needed, and activate if it is. This will restart bpp with the new control file
                # numeric scaling is not active by default because it slows bpp considerably.
                if "[ERROR] log-L for locus" in output_line:
                    print('Restarting BPP with numerical scaling')
                    file1 = open(control_file, "a")  # append mode
                    file1.write("scaling=1")
                    file1.close()
                    restart = True
                    break

                # check if bpp errored with a given error message
                elif "[ERROR]" in output_line:
                    print("#", output_line)
                    sys.exit(f"BppError: BPP failed. Check control file independently using the 'bpp' command")
                
                # check if bpp gave a segfault
                elif "core dumped" in output_line:
                    print("#", output_line)
                    sys.exit("BppError: BPP failed with segfault. Check control file independently using the 'bpp' command, and contact BPP developers if issue persists")

                # print the current progress indicator
                progress = re.findall("-*\d\d%", output_line)
                if len(progress) == 1:
                    print(f'BPP progress: {progress[0]}        ', end='\r')

            # exit if process has stopped
            if (output_line == '') and (process.poll() != None):
                break
        
        if restart == False:
            process.wait() # check again that process has stopped
            print("> Finished BPP run                             ")
            bpp_completed = True