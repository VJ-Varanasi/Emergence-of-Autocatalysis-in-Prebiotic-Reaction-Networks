
# Estimates Log of summation
def log_sum (vals):
    if vals:
        x = max(vals)
        adj_vals = [i- x for i in vals]
        log_s = x + np.log(np.sum(np.exp(adj_vals)))
        return log_s
    else:
        return 0

# Estimates size of |B_k| (eq 11)
def raf_term(k,mR, Q, k_prime = 163, n= 6, A1= 90 ):
    C = Q * np.log(k_prime) - 2* n * k_prime *np.log(2)
    #bks = 0
    bks = []

    for m in range(1, k): 
        part1 = rest_part(Q,m, C)
        singles = wiki_approx(mR - A1 - m, k- m)

        
        correction = rest_part(Q,k-m, C)

        if correction < singles:

            if (correction - singles) > -3:
                print("Not Simple")
                bks.append (np.log((np.exp(part1) * (np.exp(singles) - np.exp(correction)))))

            else:
                bks.append (part1 + singles)

        else:
            bks.append(0)
    
    if bks:

        return bks
    else:
        return 

# Combinatorial Approximation
def wiki_approx(n,k):
    if k == 0:
        return 0
    return k * (np.log(n) +1  - np.log(k))


# Restricted Parition Function Estimate (eq 14)
def rest_part(Q, n,C):
    return Q * np.log(n) - C

