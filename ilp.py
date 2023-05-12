from math import floor
import numpy as np
from fractions import Fraction as Frac

Frac.__repr__ = Frac.__str__

def swap_basic(T, i, l):
    m,n = T.shape
    n-=1
    m-=1
    T[i+1, :] = T[i+1, :]/T[i+1, l+1]
    for r in range(0, i+1):
        T[r, :] = T[r, :] - (T[r, l+1]*T[i+1, :])
    for r in range(i+2, m+1):
        T[r, :] = T[r, :] - (T[r, l+1]*T[i+1, :])
        
def simplex_tableau(T, basic):
    '''Solves the standard form Linear Programming Problem with A as the coefficient matrix
        Given: c, A, b and initial BFS x_0
        Assumes A is full rank
    '''
    #T represents some tableau and basic are the indices of the variables that are basic initially (0 indexed)
    rc = T[0, 1:] #reduced costs
    m,n = T.shape
    n-=1
    m-=1
    while(not (rc >=0).all()):
        l = np.argmin(rc>=0) #incoming variable
        '''Follows Bland's rule as we take the lowest index with reduced cost negative'''
        xB = T[1:, 0]
        dist = np.Inf
        i = None
        for r in range(m):
            if T[r+1, l+1] > 0 and xB[r]/T[r+1, l+1] < dist:
                i = r
                dist = xB[r]/T[r+1, l+1]

        if dist is np.Inf:
            #unbounded
            pass
        else:
            swap_basic(T, i, l)
            basic[i] = l

    return T, basic

def dual_simplex(T, basic):
    '''Solves the standard form Linear Programming Problem with A as the coefficient matrix
        Given: c, A, b and initial BFS x_0
        Assumes A is full rank
    '''
    #T represents some tableau and basic are the indices of the variables that are basic initially (0 indexed)

    xB = T[1:,0] #reduced costs
    m,n = T.shape
    n-=1
    m-=1
    while(not (xB >= 0).all()):
        i = np.argmin(xB >= 0) #outgoing variable
        '''Follows Bland's rule as we take the lowest index with reduced cost negative'''

        rc = T[0, 1:]
        dist = np.Inf
        l = None
        for c in range(n):
            if T[i+1][c+1] < 0 and (-rc[c]/T[i+1, c+1]) < dist:
                # print(-rc[c]/T[i+1, c+1], "Whooo")
                dist = -rc[c]/T[i+1, c+1]
                l = c

        if dist is np.Inf:
            #unbounded
            pass
        else:
            #Change the tableau
            swap_basic(T, i, l)
            basic[i] = l
            xB = T[1:,0]
    return T, basic

def isInt(a):
    return a.denominator == 1

def fractional(a):
    return a - floor(a)

def gomory(filename):
    with open(filename) as f:
        lines = f.readlines()
        n,m = map(int, lines[0].split())
        b = list(map(int, lines[1].split()))
        c = list(map(int, lines[2].split()))

        A = []
        for rr in range(3, m+3):
            A.append(list(map(int, lines[rr].split())))

        b = np.transpose(np.array(b)).reshape((-1, 1)) + Frac()
        c = np.transpose(np.array(c)).reshape((-1, 1)) + Frac()
        c = -c
        A = np.array(A) + Frac()
        
    #Add Slack vars
        ''' Variables with indices 0 to n-1 are the main variables, n to n+m-1 are the slack variables'''
        A1 = np.concatenate((A, np.identity(m, dtype=Frac)), axis=1)
        c2 = np.concatenate((c, np.zeros((m,1), dtype=Frac)+Frac()), axis=0)
    #Finding initial BFS using 2 PHASE SIMPLEX
        '''Variables with indices n+m to n+2m-1 are artificial'''
        temp = np.concatenate((A1, b), axis=1) + Frac()

        for i in range(m):
            if(b[i] < 0):
                temp[i, :] = -temp[i, :]

        A2 = np.concatenate((temp[:, 0:n+m], np.identity(m, dtype=Frac)), axis=1) + Frac()

        part1 = np.concatenate((np.array([[0]]),temp[:,-1].reshape(-1, 1)), axis = 0)
        
        newc = np.concatenate((np.zeros((n+m,1), dtype=Frac) + Frac(), np.ones((m, 1), dtype=Frac)+Frac()))
        rc1 = np.zeros((n+m, 1), dtype=Frac) + Frac()
        for i in range(n+m):
            rc1[i][0] -= np.sum(A2[:, i])
        newrc = np.concatenate((rc1, np.zeros((m,1), dtype=Frac))).reshape(-1, 1) + Frac()

        part2 = np.concatenate((np.transpose(newrc), A2), axis=0)
        T = np.concatenate((part1, part2), axis=1)
        T = T + Frac()
        T[0][0] = -np.sum(temp[:,-1].reshape(-1, 1))
        T, basic = simplex_tableau(T, np.arange(n+m, n+2*m))
        
        #Swap artificials with non-artificials
        while((basic >= n+m).any()):
            i = np.argmax(basic >= n+m) #the least index artificial variable which is in basis
            non_arti = np.arange(n+m)
            l = np.argmin(np.array([tt in basic for tt in non_arti]))
            swap_basic(T, i, l)
            basic[i] = l
        x = np.zeros((n+2*m, 1), dtype=Frac).reshape((-1, 1)) + Frac()
        x[basic, 0] = np.transpose(T[1:, 0])
        
    #Construct initial tableau
        T = T[:, 0:n+m+1]
        x = x[0:n+m, 0]

        rc = np.transpose(c2) - np.matmul(np.transpose(c2[basic,]), T[1:, 1:])
        T[0, 1:] = rc
        T[0][0] = -np.dot(x, c2)[0]

        T, basic = simplex_tableau(T, basic)

        done = 0
        while(not done):
            xB = T[1:, 0]
            i = None
            m1, n1 = T.shape
            n1 -= 1
            m1 -= 1
            for r in range(m1):
                if not isInt(xB[r]):
                    i = r
                    break
            if i is None:
                done = 1
                break
            row = T[i+1, :]
            constraint = -np.array((row%1)).reshape(1, -1)
            T = np.concatenate((T, constraint), axis=0) + Frac()
            T = np.concatenate((T, np.zeros((m1+2, 1), dtype=Frac)), axis=1) + Frac()
            T[-1][-1] = 1

            basic = np.append(basic,np.array([n1]))
            T, basic = dual_simplex(T, basic)

        m1,n1 = T.shape
        m1-=1
        n1-=1
        x = np.zeros((n1,1), dtype=Frac) + Frac()
        x[basic, 0] = T[1:, 0]
        sol = [0 for _ in range(n)]
        for k in range(n):
            sol[k] = int(x[k][0])
        return sol
    
if __name__ == "__main__":
    print(gomory("input.txt"))
