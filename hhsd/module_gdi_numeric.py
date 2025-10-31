'''
FUNCTIONS FOR CALCULATING THE GDI NUMERICALLY
'''

from .customtypehints import NodeName
from .module_ete3 import TreeNode
from .module_bpp_readres import MSCNumericParamEstimates, NumericParam

import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad


def pg1a_numeric_formula(
        theta_A:    float,
        theta_B:    float,
        tau_AB:     float,
        wAB:        float,
        wBA:        float
        ) ->        float:

    '''
    Numerically calculate P(G1) for node pairs with only reciprocal migration or no migration.
    '''

    cA = 2/theta_A
    cB = 2/theta_B
    #wAB = (4*M_AB)/theta_B   # often used entries
    #wBA = (4*M_BA)/theta_A
    

    # some rate matrix diagonals that are set up to make the rows sum to 0
    r01 = -3*(cA + wBA)       
    r02 = -(cA + wAB + 2*wBA) 
    r03 = -(cB + 2*wAB + wBA) 
    r04 = -3*(cB + wAB)       
    r05 = -2*wBA -  cA        
    r06 = -2*wAB - cB         
    r07 = - wAB - wBA        


    # rate matrix
    rate_matrix = [
    [ r01, wBA, wBA,   0, wBA,   0,   0,   0,  cA,  cA,  cA,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [ wAB, r02,   0, wBA,   0, wBA,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  cA,   0],
    [ wAB,   0, r02, wBA,   0,   0, wBA,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  cA,   0,   0],
    [   0, wAB, wAB, r03,   0,   0,   0, wBA,   0,   0,   0,   0,   0,   0,  cB,   0,   0,   0,   0,   0,   0],
    [ wAB,   0,   0,   0, r02, wBA, wBA,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  cA,   0,   0,   0],
    [   0, wAB,   0,   0, wAB, r03,   0, wBA,   0,   0,   0,   0,   0,   0,   0,  cB,   0,   0,   0,   0,   0],
    [   0,   0, wAB,   0, wAB,   0, r03, wBA,   0,   0,   0,   0,   0,   0,   0,   0,  cB,   0,   0,   0,   0],
    [   0,   0,   0, wAB,   0, wAB, wAB, r04,   0,   0,   0,  cB,  cB,  cB,   0,   0,   0,   0,   0,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0, r05,   0,   0,   0,   0,   0, wBA,   0,   0, wBA,   0,   0,  cA],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0, r05,   0,   0,   0,   0,   0, wBA,   0,   0, wBA,   0,  cA],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, r05,   0,   0,   0,   0,   0, wBA,   0,   0,  wBA, cA],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, r06,   0,   0, wAB,   0,   0, wAB,   0,   0,  cB],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, r06,   0,   0, wAB,   0,   0,  wAB,  0,  cB],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, r06,   0,   0, wAB,   0,   0, wAB,  cB],
    [   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0, r07,   0,   0,   0,   0,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0, r07,   0,   0,   0,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0, r07,   0,   0,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0,   0,   0,   0, r07,   0,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0,   0,   0,   0, r07,   0,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, wAB,   0,   0, wBA,   0,   0,   0,   0,   0, r07,   0],
    [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0]
    ]
    
    rate_matrix = np.asarray(rate_matrix, dtype=np.float64)

    
    # Get probability density of a1-a2 coalescence before tau_AB [also known as P(G1A)]
    def coalesce_prob_dens_at_tau(tau):
        
        # exponentiate matrix
        prob_matrix = expm((rate_matrix * tau))
        
        # get probability of a1-a2 coalescence events
        prob_coalescence_in_pop_A = (prob_matrix[1, 0] + prob_matrix[1, 1])*(2/theta_A)
        prob_coalescence_in_pop_B = (prob_matrix[1, 6] + prob_matrix[1, 7])*(2/theta_B)
        
        return prob_coalescence_in_pop_A + prob_coalescence_in_pop_B

    # integrate from 0 to tau_AB to get final probability of a1-a2 coalescence before tau_AB
    prob_coal_before_tau, err = quad(coalesce_prob_dens_at_tau, 0, tau_AB)
    
    pg1a = prob_coal_before_tau

    return pg1a


def get_pg1a_numerical(
        node:           TreeNode,          
        numeric_param:  MSCNumericParamEstimates,
        ) ->            NumericParam:
    
    '''
    Get the gdi of a given leaf node in the Tree object, if it can be calulcated analytically.
    Perform the 1000 replicate similations needed to establish a distribution of gdi values
    '''

    main_node:NodeName     = str(node.name)
    sister_node:NodeName   = str(node.get_sisters()[0].name)
    ancestor_node:NodeName = str(node.up.name)
    
    # perform the replicate simulations
    results = []
    for i in range(1000):
        
        print(f"inferring gdi for '{node.name}' using analytical formula ({i+1}/1000)...                    ", end = '\r')

        tau_dict        = numeric_param.sample_tau(i)
        theta_dict      = numeric_param.sample_theta(i)
        migration_df    = numeric_param.sample_migparam(i)

        # get input values 
        theta_A = theta_dict[main_node]
        theta_B = theta_dict[sister_node]
        tau_AB  = tau_dict[ancestor_node]
        
        try: # try accepts are for cases where one or both populations do not have migration to the other
            w_AB = migration_df[migration_df["source"] == main_node]['W'].to_list()[0]
        except:
            w_AB = 0

        try:
            w_BA = migration_df[migration_df["source"] == sister_node]['W'].to_list()[0]
        except:
            w_BA = 0

        # perform the calculations
        results.append(pg1a_numeric_formula(theta_A, theta_B, tau_AB, w_AB, w_BA))

    print("                                                                                                  ", end='\r')


    return NumericParam(results)