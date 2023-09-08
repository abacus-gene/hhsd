'''
FUNCTIONS FOR CALCULATING THE GDI NUMERICALLY
'''

import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad

# numerically calculate the gdi for pairs with reciprocal or no migration.
def gdi_numeric(
        theta_A,
        theta_B,
        tau_AB,
        M_AB,
        M_BA
        ):

    wAB = (4*M_AB)/theta_B   # often used entries
    wBA = (4*M_BA)/theta_A
    cA = 2/theta_A
    cB = 2/theta_B

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

    
    # probability density of a1-a2 coalescence before tau_AB
        # from eq. 10 
        # = [p_AAB,AAA(t)+p_AAB,AAB(t)]*[2/theta_A] + [p_AAB,AAA(t)+p_AAB,AAB(t)]*[2/theta_A]
    def coalesce_prob_dens_at_tau(tau):
        
        # exponentiate matrix
        prob_matrix = expm((rate_matrix * tau))
        
        # get probability of a1-a2 coalescence events
        prob_coalescence_in_pop_A = (prob_matrix[1, 0] + prob_matrix[1, 1])*(2/theta_A)
        prob_coalescence_in_pop_B = (prob_matrix[1, 6] + prob_matrix[1, 7])*(2/theta_B)
        
        return prob_coalescence_in_pop_A + prob_coalescence_in_pop_B


    '''
    gdi is the probability that the first coalescence is between the two A sequences and it 
    happens before reaching species divergence when we trace the genealogy backwards in time.
    '''
    # integrate from 0 to tau_AB to get final probability of a1-a2 coalescence before tau_AB
    res, err = quad(coalesce_prob_dens_at_tau, 0, tau_AB)

    return np.round(res, 2)

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


# wrapper function for getting adapting the numeric gdi function to the tree datastructure
def get_gdi_numerical(
        node,
        migration_df
        ):

    sister_node   = node.get_sisters()[0]
    ancestor_node = node.up

    # try accepts are for cases where one or both populations do not have migration to the other
    try:
        Mig_AB = migration_df[migration_df["source"] == str(node.name)]["M"].to_list()[0]
    except:
        Mig_AB = 0

    try:
        Mig_BA = migration_df[migration_df["source"] == str(sister_node.name)]["M"].to_list()[0]
    except:
        Mig_BA = 0

    return gdi_numeric(node.theta, sister_node.theta, ancestor_node.tau, Mig_AB, Mig_BA)