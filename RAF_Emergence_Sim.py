import numpy as np
import pandas as pd
import os, itertools, string, multiprocessing as mp
from tqdm import tqdm



# Repeat of functions presented in Kauffman_network_helper.py
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
                    #print(list(R.values()))
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

# Repeatedly applies reduction rule R1 until no more reductions can be made
def ReduceToRa(R,C):
    change = True
    while change:
        change = False
        supp = supportR(R)
        for r in R:
            # see if r is catalyzed by anything in support R
            if r in C:
                if set(C[r]).intersection(supp):
                    pass
                else:
                    #print(r)
                    #print(R)
                    del R[r]
                    change = True
                    break
            else:
                #print(r)
                #print(R)
                del R[r]
                change = True
                break



    return

def supportR(R):
    supp = []
    for r in R:
        supp += R[r][0] + R[r][1]

    return set(supp)

def create_catalysts(X, react_count, f):
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

def computeClosure(F,R):
    W =set(F)
    change = True
    while change:
        change = False 
        for r in R:
            if (set(R[r][0]).issubset(W)) and (not set(R[r][1]).issubset(W)):
                W.update(R[r][1])
                change = True
    return list(W)

def ReduceToF(R, W):

    for r in list(R.keys()):
        if not set(R[r][0]).issubset(set(W)):
            del R[r]
    return 

def RAF(F,R,C):
    
    R_old = list(R.keys())
    #print("Original: {}".format(R_old))
    change = True
    while change:
        #print("R: {}".format(list(R.keys())))
        #print("ReduceToRA")
        ReduceToRa(R,C)
        #print("ComputeClosure")
        W = computeClosure(F,R)
        #print("ReduceToF")
        ReduceToF(R,W)

        #print("New R: {}".format(list(R.keys())))
        #how do i know when to stop
        if not R:
            #print("EMPTY")
            return 0
        elif R_old == list(R.keys()):
            #print("NO CHANGE")
            change = False
            return 1
        else:
            #print("KEEP GOING")
            R_old = list(R.keys())
            


    #print("Final")
    #print(R)
    if R:
        return 1
    else:
        return 0


# ---- Simulation parameters ----
N = 25            # number of trials per setting (can vary if needed)
ns = [10]      # network sizes
fs = np.linspace(0, 2.5, 10)  # catalysis parameters

SAVE_PATH = "raf_results.csv"

# ---- Load existing results if file exists ----
if os.path.exists(SAVE_PATH):
    results_df = pd.read_csv(SAVE_PATH)
else:
    results_df = pd.DataFrame(columns=["n", "f", "trial_count", "raf_prob"])

        
# -----------------------------
# Single simulation function
# -----------------------------
def single_run(args):
    """Perform a single RAF simulation."""
    X_copy, F_copy, R_copy, f = args
    X, F, R = X_copy.copy(), F_copy.copy(), R_copy.copy()
    C = create_catalysts(X, len(R), f)
    return RAF(F, R, C)

# -----------------------------
# Parallel sweep for each n
# -----------------------------
if __name__ == "__main__":
    # Optional: this ensures safe multiprocessing on Windows/macOS
    mp.freeze_support()

    for n in ns:
        X, F, R = create_XFR(n)
        X_copy, F_copy, R_copy = X.copy(), F.copy(), R.copy()

        for f in fs:
            print(f"\nRunning simulations for n={n}, f={f:.3f} ...")

            # Prepare argument list for parallel pool
            args_list = [(X_copy, F_copy, R_copy, f)] * N

            # Run new batch of simulations
            with mp.Pool(processes=max(mp.cpu_count() - 1, 1)) as pool:
                results = list(tqdm(pool.imap_unordered(single_run, args_list), total=N))

            raf_prob_new = np.mean(results)

            # --- Check if this (n,f) pair already exists ---
            existing_mask = (results_df["n"] == n) & (np.isclose(results_df["f"], f, atol=1e-8))
            
            if existing_mask.any():
                # Retrieve existing data
                existing_row = results_df.loc[existing_mask].iloc[0]
                old_prob = existing_row["raf_prob"]
                old_N = existing_row["trial_count"]

                # Compute weighted running average
                total_N = old_N + N
                combined_prob = (old_prob * old_N + raf_prob_new * N) / total_N

                # Update row
                results_df.loc[existing_mask, "raf_prob"] = combined_prob
                results_df.loc[existing_mask, "trial_count"] = total_N

                print(f"Updated existing (n={n}, f={f:.3f}): RAF_prob={combined_prob:.3f} "
                    f"(+{N} runs, total={total_N})")

            else:
                # Add new row
                new_row = {"n": n, "f": f, "trial_count": N, "raf_prob": raf_prob_new}
                results_df = pd.concat([results_df, pd.DataFrame([new_row])], ignore_index=True)
                print(f"Added new entry: n={n}, f={f:.3f}, RAF_prob={raf_prob_new:.3f}")

            # Save incrementally
            results_df.to_csv(SAVE_PATH, index=False)


    print("Simulation complete. Results saved to:", SAVE_PATH)

