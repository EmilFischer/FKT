import numpy as np
import networkx as nx
from collections import deque

def sparsify(A, sparse) -> np.array:
    #Create dict of the form: row index -> column indices with non-zero entry.
    #Example: F[42] -> [3, 8, 14], meaning row 42 has a non-zero at column 3, 8, and 14.
    F = dict()
    lastRow = 0
    rows = deque()

    R = sparse.tocoo()
    n = A.shape[0]-1
    order = n+1

    for (r, c) in zip(R.row, R.col):
        if r not in F:
            F[r] = deque()
        F[r].append(c)

        if lastRow != r:
            if len(F[lastRow]) > 3:
                rows.append(lastRow)
                order += (len(F[lastRow])-3)*2
        lastRow = r
        
    if len(F[lastRow]) > 3:
        rows.append(r)
        order += (len(F[lastRow])-3)*2

    #Pad B with zeros for sparsification
    B = np.zeros((order, order))
    B[:-order+n+1, :-order+n+1] = A

    #Begin actual sparsification
    while len(rows) > 0:
        i = rows[0]

        #Identify two destinct indices that contain a non-zero in row i
        u = F[i].popleft()
        v = F[i].popleft()

        #Set to zero and hold old values
        B_iu = B[i, u]
        B_iv = B[i, v]
        B[i, u] = 0
        B[i, v] = 0

        B_ui = B[u, i]
        B_vi = B[v, i]
        B[u, i] = 0
        B[v, i] = 0

        #Set n+1 and n+2 entries
        B[n+1, i]   = -1
        B[n+1, n+2] = 1
        
        B[n+2, u]   = B_iu
        B[n+2, v]   = B_iv
        B[n+2, n+1] = -1

        B[i, n+1] = 1
        
        B[u, n+2] = B_ui
        B[v, n+2] = B_vi

        #Update rows i, u, and v in F
        F[i].append(n+1)
        
        if u in F:
            F[u].remove(i)
            F[u].append(n+2)

        if v in F:
            F[v].remove(i)
            F[v].append(n+2)

        #Go to next row if it contains < 4 non-zeros
        if len(F[i]) < 4:
            rows.popleft()

        #Set the new size of the matrix for next iteration
        n += 2
    
    return B

#G = nx.grid_2d_graph(100, 100)
#M = nx.to_scipy_sparse_matrix(G)
#B = sparsify(M.todense(), M)

#n = B.shape[0]
#print("n:", n)

#For checking all rows/columns have <= 3 non-zeros
#for i in range(n):
#    c = 0
#    for j in range(n):
#        if B[i, j] != 0:
#            c += 1
#    if c > 3:
#        print(i, c)