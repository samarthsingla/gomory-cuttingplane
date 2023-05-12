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


def print_tableau(T):
    for i in range(len(T)):
        for j in range(1):
            print(f"{T[i][j].numerator}/{T[i][j].denominator}", end=" | ")       
        for j in range(1,len(T[i])):
            print(f"{T[i][j].numerator}/{T[i][j].denominator}", end=" ")
        print()

# if __name__ == "__main__":
#     T = np.array([
#         [0, -10, -12, -12, 0,0,0],
#         [20, 1, 2, 2, 1, 0, 0],
#         [20, 2, 1, 2, 0, 1, 0],
#         [20, 2, 2, 1, 0, 0, 1]
#     ])
#     T = T + Frac()
#     T2, basic = simplex_tableau(T, np.array([3,4,5]))
#     print_tableau(T2)
#     print(basic)