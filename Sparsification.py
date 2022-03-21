import numpy as np
import networkx as nx
from collections import deque
import time

def sparsify(B, sparse) -> np.array:
    #Create dict of the form: row index -> column indices with non-zero entry.
    #Example: F[42] -> [3, 8, 14], meaning row 42 has a non-zero at column 3, 8, and 14.
    F = dict()
    lastRow = -1
    columns = deque()
    rowSet = set()
    rowQ = deque()
    R = sparse.tocoo()
    n = np.shape(B)[0]-1
    BOrder = n

    for (r, c) in zip(R.row, R.col):
        if lastRow != r:
            if len(columns) > 3:
                F[lastRow] = columns
                rowSet.add(lastRow)
                rowQ.append(lastRow)
                BOrder += (2*len(columns))-3
            columns = deque()
            lastRow = r

        columns.append(c)

    if len(columns) > 3:
        F[lastRow] = columns
        rowSet.add(lastRow)
        rowQ.append(lastRow)
        BOrder += len(columns)-3

    zeros = np.zeros((BOrder, BOrder))
    zeros[:-BOrder+n, :-BOrder+n] = B
    B = zeros
    
    #Begin actual sparsification algorithm
    lastRow = -1
    nPlus1 = deque()
    nPlus2 = deque()

    while len(rowSet) > 0:
        i = rowQ[0]

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
        
        t0 = time.time()
        #Add two rows and columns (n = size of B-1, cus of indexing)
        zeros = np.zeros(n+1)
        B = np.c_[B, zeros, zeros] #Adding two zero columns

        zeros = np.zeros(n+3)
        B = np.append(B, [zeros, zeros], axis=0) #Adding two zero rows
        t1 = time.time()
        total = t1-t0
        if total >= 1:
            print("Adding rows and cols:", total)

        F[n+1] = deque()
        F[n+2] = deque()

        #Set newly added entries
        B[n+1, i]   = -1
        B[n+1, n+2] = 1
        
        B[n+2, u]   = B_iu
        B[n+2, v]   = B_iv
        B[n+2, n+1] = -1

        B[i, n+1] = 1
        
        B[u, n+2] = B_ui
        B[v, n+2] = B_vi

        #Update F according to changes in row i, u, v, n+1, and n+2
        F[i].append(n+1)
        
        if u in rowSet:
            F[u].remove(i)
            F[u].append(n+2)

        if v in rowSet:
            F[v].remove(i)
            F[v].append(n+2)

        #Update F with rows n+1 and n+2 if it contains > 3 non-zeros
        if lastRow != i:
            if len(nPlus1) > 3:
                F[n+1] = nPlus1
                rowQ.append(n+1)
                rowSet.add(n+1)

                F[n+2] = nPlus2
                rowQ.append(n+2)
                rowSet.add(n+2)

            nPlus1 = deque()
            nPlus2 = deque()
            lastRow = i

        #Hold non-zero elements in rows n+1 and n+2
        nPlus1.append(i)
        nPlus1.append(n+2)

        nPlus2.append(u)
        nPlus2.append(v)
        nPlus2.append(n+1)

        #Remove # row indices with > 3 non-zeros
        if len(F[i]) < 4:
            rowSet.remove(i)
            rowQ.popleft()

        #Set the new size of the matrix for next iteration
        n += 2

    return B

G = nx.grid_2d_graph(25, 25)
M = nx.to_scipy_sparse_matrix(G)
B = sparsify(M.todense(), M)
shape = np.shape(B)
n = shape[0]
order = n + 2*M.shape[0]
print(order)
print(shape)

#for i in range(n):
#    c = 0
#    for j in range(n):
#        if B[i, j] != 0:
#            c += 1
#    if c > 3:
#        print(i, c)