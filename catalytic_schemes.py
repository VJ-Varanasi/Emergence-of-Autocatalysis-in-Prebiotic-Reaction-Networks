import numpy as np
import string
import itertools
from scipy.stats import binom

#Returns dictionaries corresponding to the molecule set, food set, and reaction set as defined by a Kauffman Network with network size n
def create_XFR(n, t = 2,k = 2):
    if n ==2:
        t = 1

    X = []
    F = []
    alphabet = string.ascii_uppercase[0:k]
    for i in range(1, n+1):    
        vals = [''.join(m) for m in itertools.product(alphabet, repeat=i)]
        X = X + vals
        if i <= t:
            F = F + vals
    
    R = {}
    react_count = 0
    for i in range(len(X)):
        cand = X[i] + X[i]
        if len(cand) <= n:
            R[react_count] = [[X[i], X[i]], [cand]]
            react_count +=1
            #Lysis Reaction
            R[react_count] = [[cand],[X[i], X[i]]]
            react_count +=1

        for j in range(i+1, len(X)):
            cand1 = X[i] + X[j]
            cand2 = X[j] + X[i]
            if len(cand1) <= n:
            
                if cand2 != cand1:
                    # Includes Reflexivity Condition
                    if [[X[j], X[i]], [cand1]] not in list(R.values()):
                        R[react_count] = [[X[i], X[j]], [cand1]]
                        react_count +=1
                    
                    if [[cand1],[X[j], X[i]]] not in list(R.values()):
                        R[react_count] = [[cand1],[X[i], X[j]]]
                        react_count +=1
                    
                    if [[X[i], X[j]], [cand2]] not in list(R.values()):
                        R[react_count] = [[X[j], X[i]], [cand2]]
                        react_count +=1
                    
                    if [[cand2],[X[i], X[j]]] not in list(R.values()):
                        R[react_count] = [[cand2],[X[j], X[i]]]
                        react_count +=1
                else:
                    if [[X[j], X[i]], [cand1]] not in list(R.values()):
                        R[react_count] = [[X[i], X[j]], [cand1]]
                        react_count +=1
                    

                    if [[cand1],[X[j], X[i]]] not in list(R.values()):
                        R[react_count] = [[cand1],[X[i], X[j]]]
                        react_count +=1
    
    return(X,F, R)

# Class used to produce theoretical estimates of Uniform Catalysis Scheme
class uniformCatalysis:
    def __init__(self, n, f):
        self.n = n 
        self.f = f
        X, F, R = create_XFR(n)  # Unpack the return values
        self.len_X = len(X)
        self.len_F = len(F)
        self.len_R = len(R)
        self.X = X  # Store the actual variables if needed
        self.F = F
        self.R = R
        self.mR = int((self.len_X * self.len_R)/2)

        ## mu  =f * X, mu = len_X * len_R/2 * p
        self.p = f/(self.len_R/2) 

        self.C_mu = self.len_X * self.len_R/2 * self.p
        self.C_var = self.len_X * self.len_R/2 * self.p * (1-self.p)

        self.mu = self.C_mu
        self.std = np.sqrt(self.C_var)

        # Plotting Helpers
        self.x_min = self.C_mu - 1/4 * self.C_var
        self.x_max = self.C_mu + 1/4 * self.C_var

    
    def cdf(self, k):
        return binom.cdf(k, self.mR, self.p)

    def pdf(self, k):
        return binom.pmf(k, self.mR, self.p)



# Class used to produce theoretical estimates of PowerLaw Catalysis Scheme (GPT 5 Assistance Used)
class powerlawCatalysis:
    def __init__(self, n, f, N=10000):
        self.n = n
        self.f = f

        X, F, R = create_XFR(n)  # assumes defined externally
        self.len_X = len(X)
        self.len_F = len(F)
        self.len_R = len(R)
        self.X, self.F, self.R = X, F, R
        self.mR = int((self.len_X * self.len_R) / 2)

        self.s = self.zipf_s_from_mean_truncated(2*f, int(self.len_R / 2))

        # mean/std for sum of shifted variables
        self.mu = f * self.len_X
        self.std = self.truncated_zipf_std(self.s, int(self.len_R / 2), self.len_X)
        self.x_min = 0 #self.mu - 2 * self.std
        self.x_max = self.mu + 1 * self.std

    # --------------------------------------------------------------------- #
    def _H_K_a(self, K, a):
        if a == 0:
            return float(K)
        k = np.arange(1, K + 1, dtype=np.float64)
        return np.sum(k ** (-a))

    def _zipf_mean_trunc(self, K, s):
        return self._H_K_a(K, s - 1) / self._H_K_a(K, s)

    def zipf_s_from_mean_truncated(self, mu_target, K, tol=1e-10, max_iter=200):
        if K < 1:
            raise ValueError("K must be a positive integer.")

        mu_target = max(1.01, mu_target)
        mu_max = K / self._H_K_a(K, 1.0)
        if not (1.0 < mu_target < mu_max):
            raise ValueError(f"Target mean must be in (1, {mu_max:.12g}) for K={K}, got {mu_target}.")

        lo, hi = 1.0 + 1e-12, 50.0
        while self._zipf_mean_trunc(K, hi) > mu_target:
            hi *= 2
            if hi > 1e6:
                break

        mu_lo = self._zipf_mean_trunc(K, lo)
        mu_hi = self._zipf_mean_trunc(K, hi)
        if not (mu_hi < mu_target < mu_lo):
            lo = 1.0 + 1e-9
            mu_lo = self._zipf_mean_trunc(K, lo)
            if not (mu_hi < mu_target < mu_lo):
                raise RuntimeError("Failed to bracket the solution; check inputs.")

        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            mu_mid = self._zipf_mean_trunc(K, mid)
            if abs(mu_mid - mu_target) <= tol:
                return mid
            if mu_mid > mu_target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def truncated_zipf_std(self, s, K, n=1):
        k = np.arange(1, K + 1, dtype=np.float64)
        Hs = np.sum(k ** (-s))
        Hs1 = np.sum(k ** (-(s - 1)))
        Hs2 = np.sum(k ** (-(s - 2)))
        mean = Hs1 / Hs
        var = Hs2 / Hs - mean**2
        var = max(var, 0.0)
        return np.sqrt(n * var)

    # --------------------------------------------------------------------- #
    def truncated_zipf_pmf_shifted(self, s: float, K: int) -> np.ndarray:
        """
        pmf p[k] for shifted truncated Zipf: values 0..K-1 (each X-1).
        """
        if s <= 1 or K < 1:
            raise ValueError("Require s>1 and K>=1 for truncated Zipf.")
        k = np.arange(1, K + 1, dtype=np.float64)
        Z = np.sum(k ** (-s))
        p = np.zeros(K, dtype=np.float64)
        p[:] = k ** (-s) / Z   # probabilities for X=1..K
        # shift: Y = X - 1 ⇒ Y in {0,..,K-1}
        return p

    def sum_of_shifted_zipf_pmf_cdf(self, s: float, K: int, n: int):
        """
        Returns (support, pmf, cdf) for S = sum_i (X_i - 1)
        where X_i ~ truncated Zipf(s) on {1..K}.
        So support = 0..n*(K-1).
        """
        if n < 1:
            raise ValueError("n must be positive.")
        base = self.truncated_zipf_pmf_shifted(s, K)

        out_len = n * (K - 1) + 1
        N = 1 << (out_len - 1).bit_length()

        F = np.fft.rfft(base, n=N)
        F_n = F ** n
        pmf_full = np.fft.irfft(F_n, n=N)

        pmf = pmf_full[:out_len].copy()
        pmf = np.maximum(pmf, 0.0)
        pmf /= pmf.sum()

        support = np.arange(0, n * (K - 1) + 1, dtype=int)
        cdf = np.cumsum(pmf)
        return support, pmf, cdf

    # convenience wrappers
    def pdf(self, k):
        support, pmf, _ = self.sum_of_shifted_zipf_pmf_cdf(
            self.s, int(self.len_R / 2), self.len_X
        )
        if 0 <= k < len(pmf):
            return pmf[k]
        return 0.0

    def cdf(self, k):
        support, _, cdf = self.sum_of_shifted_zipf_pmf_cdf(
            self.s, int(self.len_R / 2), self.len_X
        )
        if k < 0:
            return 0.0
        if k >= len(cdf):
            return 1.0
        return cdf[k]

# Class used to produce theoretical estimates of AllOrNoneCatalysis
class AllOrNoneCatalysis:

    def __init__(self, n, f):
        self.n = n 
        self.f = f
        X, F, R = create_XFR(n)  # Unpack the return values
        self.len_X = len(X)
        self.len_F = len(F)
        self.len_R = len(R)
        self.X = X  # Store the actual variables if needed
        self.F = F
        self.R = R
        self.mR = int((self.len_X * self.len_R)/2)
        self.p = (f) / (self.len_R/2)  

        self.x_min = self.p * self.mR - 1/4 * self.mR*self.p*(1-self.p)
        self.x_max = self.p * self.mR + 1/4 * self.mR*self.p*(1-self.p)

        #self.cv =1- self.p
        self.mu = self.p * self.len_X * self.len_R/2 #* 1/self.len_R *2
        self.std = np.sqrt(self.len_X * self.p * (1 - self.p) * (self.len_R/2)**2)# * (1/self.len_R)**2)

    
    def cdf(self,k):
        return binom.cdf(k/(self.len_R/2), self.len_X, self.p)
    
    def pdf(self, k):
        if k % (self.len_R/2) != 0:
            return 0
        else:
            return binom.pmf(k/(self.len_R/2), self.len_X, self.p)


# Used to produce theoretical estimates for Sparse Catalysis
class SparseCatalysis:

    def __init__(self, n, f):
        self.n = n 
        self.f = f
        X, F, R = create_XFR(n)  # Unpack the return values
        self.len_X = len(X)
        self.len_F = len(F)
        self.len_R = len(R)
        self.X = X  # Store the actual variables if needed
        self.F = F
        self.R = R
        self.mR = int((self.len_X * self.len_R)/2)
        self.p = self.f *n / (self.len_R/2)

        self.mean_Z = self.len_X * self.p * (self.len_R/2) * 1/self.n
        self.var_Z = self.len_X*(self.p * self.len_R/2 * 1/self.n * (1-1/self.n) + self.p *(1-self.p)*(self.len_R/2 * 1/self.n)**2)
        #self.len_X* self.p * (self.len_R/2) * 1/self.n * (1 - 2/self.n) + self.len_X * self.p * (1 - self.p) * (2* self.len_R/ (2 *self.n))**2



        self.x_min = self.mean_Z - 3 * np.sqrt(self.var_Z)
        self.x_max = self.mean_Z + 3 * np.sqrt(self.var_Z)

        
        self.mu = self.mean_Z
        self.std = np.sqrt(self.var_Z)
    
    def cdf(self, k):
        m = np.arange(self.len_X + 1)                  # number of active molecules
        w = binom.pmf(m, self.len_X, self.p)                # weights P(K=k)
        c_total =(m * self.len_R/2).astype(int)                            # trials in the conditional binomial
        return np.sum(w * binom.cdf(k, c_total, 1/self.n))

    def pdf(self, k):
        m = np.arange(self.len_X + 1)                  # number of active molecules
        w = binom.pmf(m, self.len_X, self.p)                # weights P(K=k)
        c_total = (m * self.len_R/2).astype(int)                            # trials in the conditional binomial
        return np.sum(w * binom.pmf(k, c_total, 1/self.n))



# Functions used to produce catalyst sets for simulation

def create_uniform_catalysts (X, react_count, f):
    C = {}
    p = f/(react_count/2)
    for i in X:
        #print(X)
        for j in range(1, react_count,2):
            if np.random.random(1)[0] < p:
                k = -1
                if j in C.keys():
                    if i not in C[j]:
                        C[j].append(i)
                else:
                    C[j]= [i]

                if j+k in C.keys():
                    if i not in C[j+k]:
                        C[j+k].append(i)
                else:
                    C[j+k]= [i]
    return C    

def create_powerlaw_catalysts (X, react_count, s):
    C = {}
    for i in X:
        num_cats = np.random.zipf(s, 1)[0] - 1
        num_cats = min(num_cats, react_count)

        #print(num_cats)
        selections = np.random.choice(np.arange(1, react_count/2), size =int(num_cats))
        #print(selections)

        for j in selections:
            if j %2 ==1:
                k = -1
            else:
                k = 1
            if j in C.keys():
                if i not in C[j]:
                    C[j].append(i)
            else:
                C[j]= [i]

            if j+k in C.keys():
                if i not in C[j+k]:
                    C[j+k].append(i)
            else:
                C[j+k]= [i]
    return C  



def create_allornone_catalysts (X, react_count, f):
    C = {}
    p = f/(react_count/2)
    for i in X:
        if np.random.random(1)[0] < p:
        #print(X)
            for j in range(1, react_count,2):
                k = -1
                if j in C.keys():
                    if i not in C[j]:
                        C[j].append(i)
                else:
                    C[j]= [i]

                if j+k in C.keys():
                    if i not in C[j+k]:
                        C[j+k].append(i)
                else:
                    C[j+k]= [i]
    return C


def create_sparse_catalysts (X, react_count, f, n =6):
    C = {}
    p = n*f/(react_count/2)
    for i in X:
        if np.random.random(1)[0] < p:
            for j in range(1, react_count,2):
                if np.random.random(1)[0] < 1/n:
                    k = -1
                    if j in C.keys():
                        if i not in C[j]:
                            C[j].append(i)
                    else:
                        C[j]= [i]

                    if j+k in C.keys():
                        if i not in C[j+k]:
                            C[j+k].append(i)
                    else:
                        C[j+k]= [i]
    return C         