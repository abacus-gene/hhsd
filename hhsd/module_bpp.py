'''
FUNCTIONS REQUIRED TO INTERACT WITH BPP
'''

import subprocess
import pandas as pd
import re
import random
import sys
import time

from .module_helper import readlines, dict_merge, get_bundled_bpp_path, format_time
from .module_msa_imap import auto_prior, auto_nloci


# Custom styler function to format floats and tuples of floats with six decimal places
def format_float(value):
    if isinstance(value, float):
        return format(value, '.6f')
    elif isinstance(value, tuple):
        return tuple(format(x, '.6f') if isinstance(x, float) else (' ' if pd.isna(x) else x) for x in value)
    elif pd.isna(value):
        return ' '  # Replace Pandas NaN with empty string
    else:
        return value


# contains the list of parameters that need to be present in a BPP control file
default_BPP_cfile_dict = {
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

# initialize state of bpp control dict based on mcf parameters, and generate parameters automatically if none are supplied

def bppctl_init(
        mcf_parameters: dict
        ) -> dict:

    '''
    Runs at the start of the pipeline. Sets up BPP parameters that will be shared throughout the iterations.

    1) collects parameters of BPP that can be set through the master control file (eg 'seed', 'threads'...)

    2) generates certain parameters that are needed for BPP to run, but not supplied.
        - seed is generated as a random int
        - nloci is the number of loci present in the MSA
        - tau and theta priors are inferred from the MSA and the delimitaiton, using the method of Bruce Rannala
    '''

    # ingest parameters from the mcf
    bpp_cdict = dict_merge(default_BPP_cfile_dict, mcf_parameters)

    # generate values for parameters if not supplied
        # seed
    if bpp_cdict['seed'] == None:
        bpp_cdict['seed'] = random.randint(1, 100000)

        # nloci
    if bpp_cdict['nloci'] == None:
        bpp_cdict['nloci'] = auto_nloci(bpp_cdict['seqfile'])
        
        # tau and theta priors
    if bpp_cdict['tauprior'] == None or bpp_cdict['thetaprior'] == None:
        priors = auto_prior(imapfile=bpp_cdict['Imapfile'], seqfile=bpp_cdict['seqfile'], tree_newick=mcf_parameters['guide_tree'])
        if bpp_cdict['tauprior']   == None:
            bpp_cdict['tauprior']   = priors['tauprior']
        if bpp_cdict['thetaprior'] == None:
            bpp_cdict['thetaprior'] = priors['thetaprior']


    return bpp_cdict

# write dict representing the BPP ".ctl" file to a text file
def bppcfile_write(
        input_dict, 
        ctl_file_name
                    ):

    # remove parameters without values
    input_dict = {item:input_dict[item] for item in input_dict if input_dict[item] != None}

    # convert to pandas dataframe
    df = pd.DataFrame(list(input_dict.items()))
    
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
        control_file,  
        ):

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
        elapsed = "0:00"
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


'''
These functions are responsible for reading the bpp outfile to get the parmeters inferred during the A00 runs. for example:

 ...

  85%  0.72 0.52 0.30 0.29 0.30 0.31   0.0034 0.0097 0.0012  0.0023 0.0007 0.0005  0.0106 0.0208   1498.79727  -2862.44226  0:50
  90%  0.72 0.52 0.30 0.29 0.30 0.31   0.0033 0.0098 0.0012  0.0022 0.0007 0.0005  0.0103 0.0206   1529.98489  -2862.55024  0:52
  95%  0.72 0.52 0.30 0.29 0.30 0.31   0.0034 0.0098 0.0012  0.0023 0.0007 0.0005  0.0099 0.0199   1509.46866  -2862.66880  0:53
 100%  0.72 0.52 0.30 0.29 0.30 0.31   0.0034 0.0098 0.0012  0.0023 0.0007 0.0005  0.0099 0.0198   1499.93147  -2862.68170  0:54

0:54 spent in MCMC

          theta_1SBC	theta_2CBC	theta_3NCA	theta_4SCA	theta_5NBC	theta_6SCANBCNCACBCSBC	theta_7SCANBCNCACBC	theta_8SCANBCNCA	theta_9SCANBC	tau_6SCANBCNCACBCSBC	tau_7SCANBCNCACBC	tau_8SCANBCNCA	tau_9SCANBC	M_SCA->NBC	M_NBC->SCA	lnL
mean      0.003409  0.009805  0.001222  0.005034  0.002877  0.000996  0.001157  0.001081  0.002238  0.002274  0.000673  0.000457  0.000191  0.009877  0.019789  -2862.681700
median    0.003246  0.008834  0.001128  0.004350  0.002110  0.000847  0.001068  0.000858  0.001527  0.002205  0.000620  0.000425  0.000180  0.000045  0.000197  -2862.586000
S.D       0.001167  0.004583  0.000520  0.003003  0.002588  0.000556  0.000573  0.000731  0.002098  0.000507  0.000243  0.000152  0.000118  0.030576  0.041703  6.864305
min       0.001194  0.002144  0.000291  0.000944  0.000363  0.000219  0.000229  0.000155  0.000311  0.001022  0.000228  0.000126  0.000019  0.000000  0.000000  -2888.362000
max       0.009767  0.035791  0.004238  0.023384  0.026804  0.006377  0.010431  0.006791  0.018704  0.005070  0.001867  0.001282  0.000813  0.330063  0.555150  -2838.653000
2.5%      0.001691  0.003914  0.000510  0.001556  0.000672  0.000371  0.000422  0.000320  0.000478  0.001479  0.000338  0.000227  0.000028  0.000000  0.000000  -2876.582000
97.5%     0.006202  0.022366  0.002484  0.014242  0.010100  0.002477  0.002399  0.003099  0.007938  0.003460  0.001249  0.000833  0.000507  0.093088  0.163909  -2849.673000
2.5%HPD   0.001497  0.002337  0.000409  0.000944  0.000433  0.000252  0.000337  0.000240  0.000311  0.001362  0.000302  0.000201  0.000021  0.000000  0.000000  -2875.902000
97.5%HPD  0.005783  0.018911  0.002278  0.010647  0.007678  0.002069  0.002158  0.002529  0.006318  0.003258  0.001176  0.000775  0.000414  0.053133  0.121130  -2849.315000
ESS*      177.528058  41.805392  123.545643  47.345292  47.317375  490.965280  178.141753  125.411967  8.609130  87.816313  36.115525  37.688442  6.805319  457.352708  30.253631  251.110649
Eff*      0.035506  0.008361  0.024709  0.009469  0.009463  0.098193  0.035628  0.025082  0.001722  0.017563  0.007223  0.007538  0.001361  0.091471  0.006051  0.050222
List of nodes, taus and thetas:
Node (+1)       Tau      Theta    Label
0          0.000000   0.003409    SBC
1          0.000000   0.009805    CBC
2          0.000000   0.001222    NCA
3          0.000000   0.005034    SCA
4          0.000000   0.002877    NBC
5          0.002274   0.000996    SCANBCNCACBCSBC
6          0.000673   0.001157    SCANBCNCACBC
7          0.000457   0.001081    SCANBCNCA
8          0.000191   0.002238    SCANBC
'''


# get a dict of the numerical+string names to string names (see lines above)
def get_node_number_map(
        BPP_outfile,                    
        ):

    lines = readlines(BPP_outfile)
    relevant_index = lines.index("List of nodes, taus and thetas:") # find the line where the node labels are listed
    lines = lines[relevant_index+2:]

    '''
    map dict takes the form 1A:A or 9ABCD:ABCD
    '''
    map_dict_long = {f'{int(line.split()[0])+1}{line.split()[3]}':f'{line.split()[3]}' for line in lines} 
    map_dict_numeric = {f'{int(line.split()[0])+1}':f'{line.split()[3]}' for line in lines} 

    return map_dict_long, map_dict_numeric

# extract the mean and 95% hpds for all parameters from the outfile
def extract_param_estimate_from_outfile(
        BPP_outfile,
        ) -> dict:

    map_dict_long, map_dict_numeric = get_node_number_map(BPP_outfile)

    # read the line listing the mean of all paramater estimates, and turn this into a dict
    outlines = readlines(BPP_outfile)

    # find the index of the relevant lines
    means_line_index = [i for i in range(len(outlines)) if "mean      " in outlines[i]][0]
    param_line_index = means_line_index - 1
    hpd_bottom_line_index = [i for i in range(len(outlines)) if "2.5%HPD   " in outlines[i]][0]
    hpd_top_line_index = [i for i in range(len(outlines)) if "97.5%HPD  " in outlines[i]][0]

    # split the lines into pieces
    full_param_names    = outlines[param_line_index].split()[0:-1]
    means               = outlines[means_line_index].split()[1:-1]
    hpd_bottoms         = outlines[hpd_bottom_line_index].split()[1:-1]
    hpd_tops            = outlines[hpd_top_line_index].split()[1:-1]
    
    # format the pieces into lists
    param_types = []
    param_nodes = []
    param_means = []
    param_hpd_95 = []

    for param_name, mean, hpd_bottom, hpd_top in zip(full_param_names, means, hpd_bottoms, hpd_tops):
        full_param_name = str(param_name).split("_", maxsplit=1)

        # get the type of parameter
        param_type = full_param_name[0]
        param_types.append(param_type)

        # get the node(s) the parameter is inferred for (get statement is used to return unaltered results for M parameters)
        param_node = map_dict_long.get(full_param_name[1], full_param_name[1])
        param_node = map_dict_numeric.get(param_node, param_node)
        param_nodes.append(param_node)
        
        # get the mean
        param_means.append(float(mean))
        
        # get the 95% hpd
        param_hpd_95.append((float(hpd_bottom), float(hpd_top)))

    # format the lists into the dataframe
    estimated_param_df = pd.DataFrame({
        'type':param_types,
        'node':param_nodes,
        'mean':param_means,
        'hpd_95':param_hpd_95
    })

    return estimated_param_df


# get a dict of node names:tau values and node names:theta values from the BPP outfile
def extract_tau_theta_values(
        estimated_param_df: pd.DataFrame,                    
        ) -> tuple[dict[str, float], dict[str, float]]:

    tau_mean_dict = estimated_param_df.query("`type` == 'tau'").set_index('node')['mean'].to_dict()
    tau_hpd_dict = estimated_param_df.query("`type` == 'tau'").set_index('node')['hpd_95'].to_dict()
    theta_mean_dict = estimated_param_df.query("`type` == 'theta'").set_index('node')['mean'].to_dict()
    theta_hpd_dict = estimated_param_df.query("`type` == 'theta'").set_index('node')['hpd_95'].to_dict()

    # create dataframe
    df = pd.DataFrame.from_dict([theta_mean_dict, theta_hpd_dict, tau_mean_dict, tau_hpd_dict])
    df = df.transpose()
    df = df.rename({0: 'theta', 1: ' ', 2: 'tau', 3: '  '}, axis=1)
    
    # write to disk
    df.to_csv("estimated_tau_theta.csv")

    # format for printing, and print to screen
    print("\nEstimated tau and theta parameters:\n")
    for col in df.columns: df[col] = df[col].apply(format_float)
    print(df.to_string(index=True, max_colwidth=36, justify="start", na_rep=' ',))
    
    return tau_mean_dict, theta_mean_dict


# read the resulting migration rates from the outfile, and return them as a pandas dataframe
def extract_mig_param_to_df(
        estimated_param_df: pd.DataFrame,
        mig
        ):
    
    # only execute if migration was inferred
    if mig is not None:

        df = estimated_param_df.copy()
        df = df[df['type'] == 'M']
        
        # check if there are any migration parameters
        if len(df) == 0:
            return None

        else:
            df.drop('type', axis=1, inplace=True)

            # split 'node' into 'source' and 'destination'
            df['source'] = df['node'].apply(lambda x: x.split('->')[0])
            df['destination'] = df['node'].apply(lambda x: x.split('->')[1])
            df = df[['source', 'destination', 'mean', 'hpd_95']]
            df.rename(columns={'mean': 'M', 'hpd_95': '  '}, inplace=True)

            # write to disk
            df.to_csv("estimated_M.csv", index=False)

            # format for printing, and print to screen
            print("\nEstimated migration rates:\n")
            for col in df.columns: df[col] = df[col].apply(format_float)
            print(df.to_string(index=False, max_colwidth=36, justify="start"))

            # format by removing hpd column before returning
            df.drop('  ', axis=1, inplace=True)
            
            return df