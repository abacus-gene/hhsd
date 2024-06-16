
from copy import deepcopy
from typing import Tuple, Dict, Optional
import pandas as pd
import numpy as np

from .customtypehints import BppMCMCfile, BppOutfile, NodeName, MigrationRates, MCMCResults
from .module_helper import readlines



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
    relevant_index = lines.index("List of nodes, taus and thetas:") # find the line where the node labels are listed
    lines = lines[relevant_index+2:]

    '''
    map dict takes the form 1A:A or 9ABCD:ABCD
    '''
    map_dict_long = {f'{int(line.split()[0])+1}{line.split()[3]}':f'{line.split()[3]}' for line in lines} 
    map_dict_numeric = {f'{int(line.split()[0])+1}':f'{line.split()[3]}' for line in lines} 

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
    Extract summary statistics for all parameters from the raw mcmc chain
    '''

    # format the pieces into lists
    param_types = []
    param_nodes = []
    param_means = []
    param_hpd_025 = []
    param_hpd_975 = []

    for col_name, values in mcmc_chain.items():
        full_param_name = str(col_name).split("_", maxsplit=1)

        # get the type of parameter
        param_type = full_param_name[0]
        param_types.append(param_type)

        # get the node(s) the parameter is inferred for (get statement is used to return unaltered results for M parameters)
        param_node = map_dict_long.get(full_param_name[1], full_param_name[1])
        param_node = map_dict_numeric.get(param_node, param_node)
        param_nodes.append(param_node)
        
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
    Print this df, and write it to disk. 
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
    Get a dataframe of M values with the mean and hpd intervals. 
    Print this df, and write it to disk. 
    '''
    
    df = numeric_param.copy()
    df = df[df['type'] == 'M']
    
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
        df.rename(columns={'mean': 'M', 'hpd_025': '2.5% HPD', 'hpd_975': '97.5% HPD'}, inplace=True)

        # write to disk
        df.to_csv("estimated_M.csv", index=False)

        # format for printing, and print to screen
        print_df = deepcopy(df)
        print("\n> Estimated migration rates:\n")
        print(print_df.to_string(index=False, max_colwidth=36, justify="start"))

class MSCNumericParamDf(pd.DataFrame):
    """
    Dataframe alias holding numericparam objects for parameters of the MSC model
    """
    pass

def evenly_spaced_integers(n, m):
    """
    Used to generate the integer indices at which the mcmc chain will be sampled
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
    Extract summary statistics for all parameters from the raw mcmc chain
    '''

    # get the indeces at which the sequence will be downsampled to
    indices = evenly_spaced_integers(n=mcmc_chain.shape[0], m=n_subsample)

    # format the pieces into lists
    param_types = []
    param_nodes = []
    param_vals  = []

    for col_name, values in mcmc_chain.items():
        full_param_name = str(col_name).split("_", maxsplit=1)

        # get the type of parameter
        param_type = full_param_name[0]
        param_types.append(param_type)

        # get the node(s) the parameter is inferred for (get statement is used to return unaltered results for M parameters)
        param_node = map_dict_long.get(full_param_name[1], full_param_name[1])
        param_node = map_dict_numeric.get(param_node, param_node)
        param_nodes.append(param_node)
        
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

        # create the dataframe holding the summary stats
        self.param_summaries = extract_param_summaries(self.mcmc_df, self.map_dict_long, self.map_dict_numeric)
        # print and save summary stats to disk
        meanhpd_tau_theta(self.param_summaries)
        meanhpd_mig(self.param_summaries)

        # extract the traces into numericParam objects
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
        Sample a set of migration rate parameters from the 
        """
        
        df = self.param_traces.copy()
        df = df[df['type'] == 'M']
        
        # if no migration patterns were inferred, return None
        if len(df) == 0:
            return None

        # otherwise get the migration rates
        else:
            df.drop('type', axis=1, inplace=True)

            # split 'node' into 'source' and 'destination'
            df['source'] = df['node'].apply(lambda x: x.split('->')[0])
            df['destination'] = df['node'].apply(lambda x: x.split('->')[1])
            df['M'] = df['val'].apply(lambda x: x[index])
            df = df[['source', 'destination', 'M']]

            return df
