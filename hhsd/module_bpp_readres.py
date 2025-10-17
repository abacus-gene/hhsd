
from copy import deepcopy
import sys
from typing import Tuple, Dict, Optional
import pandas as pd
import numpy as np

from .customtypehints import BppMCMCfile, BppOutfile, NodeName, MigrationRates, MCMCResults
from .module_helper import readlines



'''
These functions are responsible for reading the bpp outfile to get the parmeters inferred during the A00 runs. for example:

 ...
  90%  0.61 0.28 0.32 0.33 1.00 0.29 0.29   0.0102 0.0091 0.0073  0.0099 0.0002 0.0001   13600.59302   -97867.97666  0:20
  95%  0.61 0.28 0.32 0.33 1.00 0.29 0.29   0.0102 0.0091 0.0073  0.0099 0.0002 0.0001   13393.17942   -97867.84522  0:21
 100%  0.61 0.28 0.32 0.33 1.00 0.29 0.28   0.0101 0.0091 0.0072  0.0099 0.0002 0.0001   13479.06416   -97866.99213  0:22

0:22 spent in MCMC

Node-Index  Node-Type  Node-Label
---------------------------------
1           Tip        D
2           Tip        A
3           Tip        B
4           Tip        C
5           Root       ABCD
6           Inner      ABC
7           Inner      BC

 param     mean     median      S.D       min       max       2.5%     97.5%    2.5%HPD   97.5%HPD     ESS*       Eff*      rho1  
----------------------------------------------------------------------------------------------------------------------------------
theta:1  0.010144  0.010126   0.000528  0.008406  0.011987  0.009139  0.011200  0.009109  0.011162  549.408569  0.109882  0.368798
theta:2  0.009050  0.008542   0.002574  0.004008  0.024486  0.005531  0.015671  0.005084  0.014322  191.511917  0.038302  0.834114
theta:3  0.007241  0.006815   0.002586  0.002228  0.031983  0.003453  0.013490  0.002846  0.012184  165.774384  0.033155  0.775707
theta:4  0.008289  0.007739   0.003333  0.002116  0.053735  0.003949  0.015993  0.003391  0.013954  310.056698  0.062011  0.769121
theta:5  0.008458  0.008445   0.000984  0.005723  0.012038  0.006645  0.010477  0.006549  0.010360   40.956773  0.008191  0.553858
theta:6  0.009788  0.009790   0.000424  0.008259  0.011289  0.008963  0.010621  0.008973  0.010624  157.649934  0.031530  0.604400
theta:7  0.013394  0.011047   0.010062  0.001870  0.167585  0.003559  0.036992  0.002154  0.028605  538.488202  0.107698  0.709658
tau:5    0.009888  0.009876   0.000380  0.008699  0.011214  0.009196  0.010668  0.009135  0.010589   16.131708  0.003226  0.834069
tau:6    0.000154  0.000151   0.000024  0.000105  0.000256  0.000114  0.000213  0.000108  0.000198   34.004855  0.006801  0.984226
tau:7    0.000105  0.000107   0.000028  0.000025  0.000198  0.000049  0.000161  0.000047  0.000156   87.619478  0.017524  0.836812

lnL      -97866.992134  -97866.341000  31.874816  -98009.035000  -97757.674000  -97932.002000  -97805.293000  -97927.665000  -97801.759000  270.384452  0.054077  0.637395

List of nodes, taus and thetas:
Node (+1)       Tau      Theta    Label
0          0.000000   0.010144       D [ D ]
1          0.000000   0.009050       A [ A ]
2          0.000000   0.007241       B [ B ]
3          0.000000   0.008289       C [ C ]
4          0.009888   0.008458    ABCD [ D A B C ]
5          0.000154   0.009788     ABC [ A B C ]
6          0.000105   0.013394      BC [ B C ]

Summarizing parameter estimates using file hhsd_job.conditional_a1b1.txt ...

 param     mean      S.D       2.5%     97.5%    2.5%HPD   97.5%HPD    Effu      Effy       c    
-------------------------------------------------------------------------------------------------
theta:1  0.010144  0.000528  0.009150  0.011218  0.009112  0.011176  0.134351  0.125048  2.240845
theta:2  0.009081  0.002610  0.005534  0.015626  0.004909  0.014279  0.040358  0.040083  1.204501
theta:3  0.007261  0.002601  0.003454  0.013452  0.002782  0.012180  0.034681  0.034522  1.152407
theta:4  0.008281  0.003249  0.003947  0.015712  0.003181  0.013947  0.067897  0.067247  1.165996
theta:5  0.008456  0.000977  0.006638  0.010438  0.006521  0.010313  0.008174  0.008144  1.816830
theta:6  0.009789  0.000428  0.008954  0.010637  0.008940  0.010623  0.033336  0.033043  1.363934
theta:7  0.013423  0.010014  0.003610  0.036821  0.002326  0.029281  0.112212  0.111038  1.104022
'''



def hpd(values:np.ndarray, mass_frac: float) :
    """
    Returns highest probability density region given by
    a set of samples.

    """
    # Get sorted list
    d = np.sort(np.copy(values))

    # Number of total samples taken
    n = len(values)
    
    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)
    
    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n-n_samples]
    
    # Pick out minimal interval
    min_int = np.argmin(int_width)
    
    # Return interval
    return [d[min_int], d[min_int+n_samples]]


class NumericParam():
    """
    class holding replicate results for a numeric param
    """
    def __init__(self, values: list[float]):
        self.values = np.array(values)
    
    def percentile(self, percentile):
        return np.percentile(self.values, percentile)
    
    def mean(self):
        return np.mean(self.values)
    
    def hpd_bound(self, mass_frac):
        return hpd(self.values, mass_frac)

    def __repr__(self) -> str:
        return f"{self.values.shape[0]} Values, Mean: {np.round(self.mean(), 6)}"
    

def read_bpp_mcmc_out(
        BPP_mcmcfile:   BppMCMCfile
        ) ->            MCMCResults:
    
    """
    read the raw mcmc results from bpp into a dataframe
    """
    lines = readlines(BPP_mcmcfile)
    if len(lines) < 2:
        sys.exit("Error: BPP mcmc output file is empty")
    # read the raw MCMC data file
    mcmc_chain = pd.read_csv(BPP_mcmcfile, delimiter='\t')
    # drop first and last columns corresponding to the Gen and lNL values, which are irrelevant
    mcmc_chain.drop(columns=[mcmc_chain.columns[0], mcmc_chain.columns[-1]], inplace=True)

    return mcmc_chain


def get_node_number_map(
        BPP_outfile:    BppOutfile,                    
        ) ->            Tuple[dict, dict]:
    
    '''
    Get a dict of the numerical+string names to string names for each node in the tree.\\
    This is used to map parameter values in the outfile to their specific nodes.
    '''

    lines = readlines(BPP_outfile)

    relevant_index = next((i for i, s in enumerate(lines) if s.startswith("Node-Index")), -1)
    if relevant_index == -1:
        sys.exit("Error: could not find node index section in BPP output file")

    lines = lines[relevant_index+2:]
    end_index = next((i for i, s in enumerate(lines) if 'param' in s), -1)
    lines = lines[:end_index-1]

    '''
    map dict takes the form 1A:A or 9ABCD:ABCD
    '''
    map_dict_long = None#{f'{int(line.split()[0])+1}{line.split()[3]}':f'{line.split()[3]}' for line in lines} 
    map_dict_numeric = None#{f'{int(line.split()[0])+1}':f'{line.split()[3]}' for line in lines} 

    return map_dict_long, map_dict_numeric

class MSCNumericParamSummary(pd.DataFrame):
    """
    Alias class for the dataframe holding summary statistics for numeric parameters for which MCMC sampling was performed by bpp.
    """
    pass

def extract_param_summaries(
        mcmc_chain:         MCMCResults,
        map_dict_long:      dict,
        map_dict_numeric:   dict,
        ) ->                MSCNumericParamSummary:            

    '''
    Extract summary statistics (mean, 2.5% HPD, 97.5% HPD) for all parameters from the raw mcmc chain
    '''

    # format the pieces into lists
    param_types = []
    param_nodes = []
    param_means = []
    param_hpd_025 = []
    param_hpd_975 = []

    for col_name, values in mcmc_chain.items():
        param_type, node_index, popname = str(col_name).split(":", maxsplit=2)

        param_types.append(param_type)
        param_nodes.append(popname)
        
        # get the actual values from the mcmc chain
        numericpar = NumericParam(np.array(values))

        # get the mean from the chain
        param_means.append(float(numericpar.mean()))

        # get the 2.5% and 97.5% hpd bounds
        hpd_bounds = numericpar.hpd_bound(0.95)
        param_hpd_025.append(float(hpd_bounds[0]))
        param_hpd_975.append(float(hpd_bounds[1]))


    # format the lists into the dataframe
    numeric_param_summary = pd.DataFrame({
        'type':param_types,
        'node':param_nodes,
        'mean': param_means,
        'hpd_025': param_hpd_025,
        'hpd_975':param_hpd_975
    })
    
    return numeric_param_summary

def meanhpd_tau_theta(
        numeric_param:  MSCNumericParamSummary                    
        ):

    '''
    Get a dataframe of tau and theta values with the mean and hpd intervals. 
    Reformat, print this df to the screen, and write it to disk. 
    '''

    tau_mean_dict       = numeric_param.query("`type` == 'tau'").set_index('node')['mean'].to_dict()
    tau_hpd_025_dict    = numeric_param.query("`type` == 'tau'").set_index('node')['hpd_025'].to_dict()
    tau_hpd_975_dict    = numeric_param.query("`type` == 'tau'").set_index('node')['hpd_975'].to_dict()

    theta_mean_dict     = numeric_param.query("`type` == 'theta'").set_index('node')['mean'].to_dict()
    theta_hpd_025_dict  = numeric_param.query("`type` == 'theta'").set_index('node')['hpd_025'].to_dict()
    theta_hpd_975_dict  = numeric_param.query("`type` == 'theta'").set_index('node')['hpd_975'].to_dict()

    # create dataframe
    df = pd.DataFrame.from_dict([theta_mean_dict, theta_hpd_025_dict, theta_hpd_975_dict, tau_mean_dict, tau_hpd_025_dict, tau_hpd_975_dict])
    df = df.transpose()
    df = df.rename({0: 'theta', 1: '2.5% HPD', 2: '97.5% HPD', 3: 'tau', 4: '2.5% HPD', 5: '97.5% HPD'}, axis=1)
    
    # write results to disk
    df.to_csv("estimated_tau_theta.csv")

    # format for printing, and print to screen
    print("\n> Estimated tau and theta parameters:\n")
    print(df.to_string(index=True, max_colwidth=36, justify="start", na_rep=' ',))

def meanhpd_mig(
        numeric_param:  MSCNumericParamSummary,    
        ):
    
    '''
    Get a dataframe of W values with the mean and hpd intervals. 
    Reformat, print this df to the screen, and write it to disk for the user. 
    '''
    
    df = numeric_param.copy()
    df = df[df['type'] == 'W']
    
    # if no migration patterns were inferred, return None
    if len(df) == 0:
        pass

    # otherwise get the migration rates
    else:
        df.drop('type', axis=1, inplace=True)

        # split 'node' into 'source' and 'destination'
        df['source'] = df['node'].apply(lambda x: x.split('->')[0])
        df['destination'] = df['node'].apply(lambda x: x.split('->')[1])
        df = df[['source', 'destination', 'mean', 'hpd_025', 'hpd_975']]
        df.rename(columns={'mean': 'W', 'hpd_025': '2.5% HPD', 'hpd_975': '97.5% HPD'}, inplace=True)

        # write to disk
        df.to_csv("estimated_W.csv", index=False)

        # format for printing, and print to screen
        print_df = deepcopy(df)
        print("\n> Estimated mutation scaled migration rates:\n")
        print(print_df.to_string(index=False, max_colwidth=36, justify="start"))

class MSCNumericParamDf(pd.DataFrame):
    """
    Dataframe alias holding numericparam objects for parameters of the MSC model
    """
    pass

def evenly_spaced_integers(n, m):
    """
    Used to generate the integer indices at which the mcmc chain will be sampled
    for example, if nsample = 10000 in the control file, an MCMC chain of length 10000 will be written to disk by bpp.
    If we want to get say 1000 samples, we take every 10th element of the full chain. 
    This function returns the indeces of the chosen chain elements
    """

    # Generate count evenly spaced numbers from n to m
    values = np.linspace(0, n, num=m, endpoint=False)
    # Convert to integers and handle duplicates
    return np.unique(values.round().astype(int))

def extract_param_traces(
        mcmc_chain:         MCMCResults,
        map_dict_long:      dict,
        map_dict_numeric:   dict,
        n_subsample:        int = 1000,          
        ) ->                MSCNumericParamDf:            

    '''
    Extract "n_subsample" evenly spaced samples from the full MCMC chain for the numeric parameters of each node
    these samples are later used to calculate distributions over the gdi
    '''

    # get the indeces at which the sequence will be downsampled to
    indices = evenly_spaced_integers(n=mcmc_chain.shape[0], m=n_subsample)

    # format the pieces into lists
    param_types = []
    param_nodes = []
    param_vals  = []

    for col_name, values in mcmc_chain.items():
        param_type, node_index, popname = str(col_name).split(":", maxsplit=2)

        param_types.append(param_type)
        param_nodes.append(popname)
        
        # get the actual values from the mcmc chain (thinned to n_subsample samples)
        values_arr = np.array(values)
        param_vals.append(values_arr[indices])

    # format the lists into the dataframe
    numeric_param_df = pd.DataFrame({
        'type':param_types,
        'node':param_nodes,
        'val':param_vals
    })
    
    return numeric_param_df



class MSCNumericParamEstimates():
    def __init__(self, BPP_outfile: BppOutfile, BPP_mcmcfile: BppMCMCfile):
        # read in the actual mcmc results
        self.mcmc_df = read_bpp_mcmc_out(BPP_mcmcfile)
        
        # read in some info that helps map between actual node names, and the names used by bpp (this is needed due to bpp shortening overly long species names in the outfile and mcmc)
        map_dict_long, map_dict_numeric = get_node_number_map(BPP_outfile)
        self.map_dict_long = map_dict_long
        self.map_dict_numeric = map_dict_numeric

        # create the dataframe holding the summary stats (mean and HPD intervals)
        self.param_summaries = extract_param_summaries(self.mcmc_df, self.map_dict_long, self.map_dict_numeric)
        # print and save summary stats to disk
        meanhpd_tau_theta(self.param_summaries)
        meanhpd_mig(self.param_summaries)

        # extract the traces (1000 evenly spaced samples from the MCMC chain) into numericParam objects
        self.param_traces = extract_param_traces(self.mcmc_df, self.map_dict_long, self.map_dict_numeric)

    def sample_tau(self, index:int) -> Dict[NodeName, float]:
        """
        Sample a tau value dictionary from the mcmc chain at the given index
        """
        tau_dict = {}
        for node_name, row in self.param_traces.query("`type` == 'tau'").set_index('node').iterrows():
            tau_dict[node_name] = row["val"][index]

        return tau_dict

    def sample_theta(self, index:int) -> Dict[NodeName, float]:
        """
        Sample a theta value dictionary from the mcmc chain at the given index
        """
        theta_dict = {}
        for node_name, row in self.param_traces.query("`type` == 'theta'").set_index('node').iterrows():
            theta_dict[node_name] = row["val"][index]

        return theta_dict

    def sample_migparam(self, index:int) -> Optional[MigrationRates]:
        """
        Sample a set of migration rate parameters from the mcmc chain
        """
        
        df = self.param_traces.copy()
        df = df[df['type'] == 'W']
        
        # if no migration patterns were inferred, return None
        if len(df) == 0:
            return None

        # otherwise get the migration rates
        else:
            df.drop('type', axis=1, inplace=True)

            # split 'node' into 'source' and 'destination'
            df['source'] = df['node'].apply(lambda x: x.split('->')[0])
            df['destination'] = df['node'].apply(lambda x: x.split('->')[1])
            df['W'] = df['val'].apply(lambda x: x[index])
            df = df[['source', 'destination', 'W']]

            return df
