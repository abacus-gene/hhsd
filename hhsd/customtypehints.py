'''
CUSTOM CLASSES FOR TYPE HINTS
'''

import pandas as pd
from typing import NewType, TypeVar, Literal, List


class CfileParam(dict):
    '''
    Alias class for hhsd control file parameters.
    '''
    pass

class Cfile(str):
    '''
    Alias class for the hhsd control file (or filepath).
    '''
    pass

class BppCfileParam(dict):
    '''
    Alias class for bpp control file parameters.
    '''
    pass

class BppCfile(str):
    '''
    Alias class for the bpp control file (or filepath).
    '''
    pass

class BppOutfile(str):
    '''
    Alias class for the bpp output file (or filepath).
    '''
    pass

class NumericParamEstimates(pd.DataFrame):
    '''
    Dataframe containing the tau, theta, and possibly M parameters for each node in the MSC model
    '''
    pass

class MigrationPattern(pd.DataFrame):
    '''
    Dataframe containing 'source' and 'destination' columns for migration events.
    '''
    pass

class MigrationRates(pd.DataFrame):
    '''
    Dataframe containing 'source', 'destination', and 'M' (migration rate) columns for migration events.
    '''
    pass

class NodeName(str):
    '''
    Name of a given node in the tree
    '''
    pass

gdi = NewType('gdi', float)

GeneTrees = TypeVar('GeneTrees', bound=List[str])

AlgoMode = Literal['merge', 'split']

Bound = Literal['lower', 'mean', 'upper']

class NewickTree(str):
    '''
    Alias class for a newick tree
    '''
    pass

class ImapPopInd(dict):
    '''
    Imap dict with populations as keys and lists of individuals in the population as values
    '''
    pass

class ImapIndPop(dict):
    '''
    Imap dict with individuals as keys and the population they belong to as values
    '''
    pass

class Filename(str):
    '''
    Alias class for general filenames
    '''
    pass