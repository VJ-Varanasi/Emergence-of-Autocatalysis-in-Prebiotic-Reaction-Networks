import string
import itertools

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
                    del R[r]
                    change = True
                    break
            else:
                del R[r]
                change = True
                break
    return

# Gets Support of Reaction SEt
def supportR(R):
    supp = []
    for r in R:
        supp += R[r][0] + R[r][1]
    return set(supp)

# Computes Closure given Food and reaction Set
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

# Reduces to food generated set
def ReduceToF(R, W):
    for r in list(R.keys()):
        if not set(R[r][0]).issubset(set(W)):
            del R[r]
    return 

#Implements RAF algorithm introduced in 'W. Hordijk and M. Steel, Detecting autocatalytic, selfsustaining sets in chemical reaction systems'
def RAF(F,R,C):
    R_old = list(R.keys())
    change = True
    while change:

        ReduceToRa(R,C)
        W = computeClosure(F,R)
        ReduceToF(R,W)

        if not R:
            return 0
        elif R_old == list(R.keys()):
            change = False
            return 1
        else:
            R_old = list(R.keys())
    if R:
        return 1
    else:
        return 0
   
