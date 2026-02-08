import numpy as np
from scipy.special import logsumexp
from scipy.special import gammaln
# Estimates Log of summation
# def log_sum (vals):
#     if vals:
#         x = max(vals)
#         adj_vals = [i- x for i in vals]
#         log_s = x + np.log(np.sum(np.exp(adj_vals)))
#         return log_s
#     else:
#         return 0

def log_sum(vals, nan_policy="omit"):
    v = np.asarray(vals, dtype=float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return -np.inf   # log(0)

    return logsumexp(v)

# Estimates size of |B_k| (eq 11)
def raf_term(k,mR, Q, k_prime = 163, n= 6 ):
    C = Q * np.log(k_prime) - 2* n * k_prime *np.log(2)

    #bks = 0
    bks = []

    for m in range(1, k): 
        part1 = rest_part(Q,m, C)
        singles = log_gamma_approx(mR - 96 - m, k- m)

        
        correction = log_sum([rest_part(Q,j, C) for j in range(2, k-m)])

        # Does the correction term dominate single correction term
        if correction < singles:
            
            # Only consider correction if magnitude is signficant
            if (correction - singles) > -4:
                print("Not Simple")
                corrected_term = np.log(np.exp(singles) - np.exp(correction))
                bks.append (part1 + corrected_term)

            else:
                bks.append (part1 + singles)

        # Otherwise contribute 0
        else:
            bks.append(-np.inf)
    
    if bks:

        return bks
    else:
        return 

# Combinatorial Approximation
def wiki_approx(n,k):
    if k == 0:
        return 0
    return k * (np.log(n) +1  - np.log(k))

def log_gamma_approx(n,k):
    return gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)

# Restricted Parition Function Estimate (eq 14)
def rest_part(Q, n,C):
    return Q * np.log(n) - C

