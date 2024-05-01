'''
FUNCTIONS FOR CALCULATING THE GDI NUMERICALLY
'''

from .customtypehints import gdi, MigrationRates, Bound
from .module_ete3 import TreeNode

import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad


def pg1a_numeric_formula(
        theta_A:    float,
        theta_B:    float,
        tau_AB:     float,
        M_AB:       float,
        M_BA:       float
        ) ->        float:

    '''
    Numerically calculate P(G1) for node pairs with only reciprocal migration or no migration.
    '''

    cA = 2/theta_A
    cB = 2/theta_B
    wAB = (4*M_AB)/theta_B   # often used entries
    wBA = (4*M_BA)/theta_A
    

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

    # # alternative calculation based on eigenvectors and eigenvalues (does not work when both M = 0)
    #
    # eigval, V, U = eig(rate_matrix, left=True, right=True)
    # eigval = np.real(eigval)
    # eigval[20] = -1
    # eigval = np.array([(np.exp(eigval[k]*tau_AB) - 1)/eigval[k] for k in range(0, 21)])
    # eigval[20] = tau_AB
    # P = U.dot(np.diag(eigval)).dot(np.linalg.inv(U)) 

    # prob_coalescence_in_pop_A = (P[1, 0] + P[1, 1])*(2/theta_A)
    # prob_coalescence_in_pop_B = (P[1, 6] + P[1, 7])*(2/theta_B)
    
    # print (prob_coalescence_in_pop_A + prob_coalescence_in_pop_B)


def get_pg1a_numerical(
        node:           TreeNode,          
        migration_df:   MigrationRates,
        bound:          Bound
        ) ->            float:
    
    '''
    Get the gdi of a given leaf node in the Tree object.

    node is a TreeNode object representing a specific population
    bound is the bound of the gdi value to be calculated, either 'lower', 'mean', or 'upper'
    migration_df is the dataframe holding the migration rates
    '''

    sister_node:TreeNode   = node.get_sisters()[0]
    ancestor_node:TreeNode = node.up
    
    # get input values based on the bound
    if bound == "lower":
        migrate_key = '97.5% HPD' # higher migration rates decrease the gdi
        theta_A = node.theta_hpd_975 # higher theta values lead to a lower gdi
        theta_B = sister_node.theta_hpd_975 
        tau_AB = ancestor_node.tau_hpd_025 # lower tau values lead to a lower gdi
    elif bound == "mean":
        migrate_key = 'M'
        theta_A = node.theta_mean
        theta_B = sister_node.theta_mean
        tau_AB = ancestor_node.tau_mean
    elif bound == "upper":
        migrate_key = '2.5% HPD' # lower migration rates increase the gdi
        theta_A = node.theta_hpd_025 # lower theta values lead to a higher gdi
        theta_B = sister_node.theta_hpd_025
        tau_AB = ancestor_node.tau_hpd_975 # higher tau values lead to a higher gdi

    # try accepts are for cases where one or both populations do not have migration to the other
    try:
        Mig_AB = migration_df[migration_df["source"] == str(node.name)][migrate_key].to_list()[0]
    except:
        Mig_AB = 0

    try:
        Mig_BA = migration_df[migration_df["source"] == str(sister_node.name)][migrate_key].to_list()[0]
    except:
        Mig_BA = 0

    return pg1a_numeric_formula(theta_A, theta_B, tau_AB, Mig_AB, Mig_BA)