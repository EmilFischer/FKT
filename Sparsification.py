import numpy as np
import networkx as nx

def sparsify(matrix, sparse) -> np.array:
    #Create list where index match row index and element is # non-zeros is row
    #Example: F[42] = 69, means row 42 contain 69 non-zero elements
    F = []
    lastRow = -1
    R = sparse.tocoo()
    rowsWithMoreThan3 = 0
    
    for (r, _) in zip(R.row, R.col):
        if r != lastRow:
            if len(F) > 0 and F[-1] > 3:
                rowsWithMoreThan3 += 1
            lastRow = r
            F.append(0)
        F[-1] += 1

    B = matrix.copy()
    n = np.shape(B)[0]-1
    while rowsWithMoreThan3 > 0:
        i = -1
        for k in range(n):
            if F[k] > 3:
                i = F[k]
                break
        #UPDATE F To CONTAIN CORRECT # NON-ZEROS
        #WHEN SETTING NEWLY ADDED ENTRIES!!!!!!
        
        
        #Identify two neighbors of vtx i
        u = -1
        v = -1
        for j in range(n):
            if B[i, j] != 0:
                if u == -1:
                    u = j
                else:
                    v = j
                    break

        #Set to zero and hold old values
        B_iu = B[i, u]
        B_iv = B[i, v]
        B[i, u] = 0
        B[i, v] = 0

        B_ui = B[u, i]
        B_vi = B[v, i]
        B[u, i] = 0
        B[v, i] = 0
        
        #Add two rows and columns (n = size of B-1, cus of indexing)
        B.resize(n+3, n+3)

        #Set newly added entries
        B[n+1, i]   = -1
        B[n+1, n+2] = 1
        
        B[n+2, u]   = B_iu
        B[n+2, v]   = B_iv
        B[n+2, n+1] = -1

        B[i, n+1] = 1
        
        B[u, n+2] = B_ui
        B[v, n+2] = B_vi

        #Update the number of elements in row i and check if < 4
        c = 0
        for j in range(n+2):
            if B[i, j] != 0:
                c += 1
        if c < 4:
            F.pop()

        #Set the new size of the matrix for next iteration
        n += 2

    return B

G = nx.grid_2d_graph(5, 5)
M = nx.to_scipy_sparse_matrix(G)
B = sparsify(M.todense(), M)
n = np.shape(B)[0]

for i in range(n):
    c = 0
    for j in range(n):
        if B[i, j] != 0:
            c += 1
    if c > 3:
        print(i, c)
